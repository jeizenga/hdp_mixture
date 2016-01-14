#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <inttypes.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "sonLib.h"
#include "ranlib.h"

#define N_IG_NUM_PARAMS 4
#define MINUS_INF -0.5 * DBL_MAX

typedef struct Factor Factor;
typedef struct DirichletProcess DirichletProcess;
typedef enum FactorType FactorType;

enum FactorType {
    BASE,
    MIDDLE,
    DATA_PT
};

struct Factor {
    FactorType factor_type;
    struct Factor* parent;
    stSet* children;
    double* factor_data;
    struct DirichletProcess* dp;
};

struct DirichletProcess {
    int id;
    struct HierarchicalDirichletProcess* hdp;
    double* gamma;
    int depth;

    struct DirichletProcess* parent;
    stList* children;
    stSet* factors;
    int num_factor_children;

    double base_factor_wt;
    double* posterior_predictive;
    double* spline_slopes;

    double cached_factor_mean;
    double cached_factor_sum_sq_dev;
    int cached_factor_size;

    bool observed;
};

struct HierarchicalDirichletProcess {
    bool finalized;
    double* data;
    int* data_pt_dp_id;
    int data_length;

    struct DirichletProcess* base_dp;
    struct DirichletProcess** dps;
    int num_dps;

    // normal-inverse gamma parameters
    double mu;
    double nu;
    double two_alpha;
    double beta;

    double* sampling_grid;
    int grid_length;
    
    int samples_taken;
    bool splines_finalized;
    
    // TODO: replace this with my new offset log gamma memo
    struct SumOfLogsMemo* log_sum_memo;
    
    int depth;
    bool sample_gamma;
    double* gamma;
    
    double* gamma_alpha;
    double* gamma_beta;
    double* w_aux_vector;
    bool* s_aux_vector;
};

void cache_base_factor_params(Factor* fctr, double mu, double nu, double two_alpha, double beta, double log_post_term) {
    if (fctr->factor_type != BASE) {
        fprintf(stderr, "Can only cache parameters for base factors.\n");
        exit(EXIT_FAILURE);
    }

    double* param_array = fctr->factor_data;
    param_array[0] = mu;
    param_array[1] = nu;
    param_array[2] = two_alpha;
    param_array[3] = beta;
    param_array[4] = log_post_term;
}

Factor* new_base_factor(HierarchicalDirichletProcess* hdp) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = BASE;

    fctr->factor_data = (double*) malloc(sizeof(double) * (N_IG_NUM_PARAMS + 1));
    cache_base_factor_params(fctr, hdp->mu, hdp->nu, hdp->two_alpha, hdp->beta, 1.0);

    fctr->parent = NULL;
    fctr->children = stSet_construct();

    DirichletProcess* base_dp = hdp->base_dp;
    fctr->dp = base_dp;
    stSet_insert(base_dp->factors, (void*) fctr);

    return fctr;
}

Factor* new_middle_factor(DirichletProcess* dp) {
    if (dp->parent == NULL) {
        fprintf(stderr, "Attempted to create middle factor in root Dirichlet process.\n");
        exit(EXIT_FAILURE);
    }

    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = MIDDLE;
    fctr->factor_data = NULL;

    // note: assigning to parent handled externally
    fctr->parent = NULL;
    fctr->children = stSet_construct();

    fctr->dp = dp;
    stSet_insert(dp->factors, (void*) fctr);
    return fctr;
}

Factor* new_data_pt_factor(HierarchicalDirichletProcess* hdp, int data_pt_idx) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = DATA_PT;
    fctr->factor_data = &(hdp->data[data_pt_idx]);

    // note: assigning to parent handled externally
    fctr->parent = NULL;
    fctr->children = NULL;

    fctr->dp = NULL;

    return fctr;
}

void destroy_factor(Factor* fctr) {
    stSet* children = fctr->children;
    if (children != NULL) {
        if (stSet_size(children) > 0) {
            fprintf(stderr, "Attempted to destroy factor that still has children.\n");
            exit(EXIT_FAILURE);
        }
        stSet_destruct(children);
    }

    Factor* parent = fctr->parent;
    if (parent != NULL) {
        stSet_remove(parent->children, (void*) fctr);
        (parent->dp->num_factor_children)--;
        if (stSet_size(parent->children) == 0) {
            destroy_factor(parent);
        }
    }

    if (fctr->factor_type == BASE) {
        free(fctr->factor_data);
    }
    
    DirichletProcess* dp = fctr->dp;
    if (dp != NULL) {
        stSet_remove(dp->factors, (void*) fctr);
    }

    free(fctr);
}

Factor* get_base_factor(Factor* fctr) {
    while (fctr->factor_type != BASE) {
        fctr = fctr->parent;
        if (fctr == NULL) {
            break;
        }
    }
    return fctr;
}

double get_factor_data_pt(Factor* fctr) {
    if (fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Attempted to access data point from non-leaf factor.\n");
        exit(EXIT_FAILURE);
    }
    return *(fctr->factor_data);
}

void get_factor_sum_internal(Factor* fctr, double* sum, int* num_data) {
    if (fctr->factor_type == DATA_PT) {
        *sum += get_factor_data_pt(fctr);
        // TODO: there should be a way to use the parent's counters instead of recounting the data pts
        (*num_data)++;
    }
    else {
        stSetIterator* child_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            get_factor_sum_internal(child_fctr, sum, num_data);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        stSet_destructIterator(child_iter);
    }
}

void get_factor_ssd_internal(Factor* fctr, double center, double* sum_sq_dev) {
    if (fctr->factor_type == DATA_PT) {
        double dev = get_factor_data_pt(fctr) - center;
        *sum_sq_dev += dev * dev;
    }
    else {
        stSetIterator* child_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            get_factor_ssd_internal(child_fctr, center, sum_sq_dev);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        stSet_destructIterator(child_iter);
    }
}

void get_factor_stats(Factor* fctr, double* mean_out, double* sum_sq_dev_out, int* num_data_out) {
    *mean_out = 0.0;
    *sum_sq_dev_out = 0.0;
    *num_data_out = 0;
    get_factor_sum_internal(fctr, mean_out, num_data_out);
    *mean_out /= (double) *num_data_out;
    get_factor_ssd_internal(fctr, *mean_out, sum_sq_dev_out);
}

void add_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    double* param_array = fctr->factor_data;

    double mu_prev = param_array[0];
    double nu_prev = param_array[1];
    double two_alpha_prev = param_array[2];
    double beta_prev = param_array[3];

    double nu_post = nu_prev + num_data;
    double mu_post = (mu_prev * nu_prev + mean * num_data) / nu_post;
    double two_alpha_post = two_alpha_prev + num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_post = beta_prev + .5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_post, two_alpha_post, beta_post,
                                                          fctr->dp->hdp->log_sum_memo);

    cache_base_factor_params(fctr, mu_post, nu_post, two_alpha_post, beta_post, log_post_term);
}

void remove_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    double* param_array = fctr->factor_data;

    double mu_post = param_array[0];
    double nu_post = param_array[1];
    double two_alpha_post = param_array[2];
    double beta_post = param_array[3];

    double nu_prev = nu_post - num_data;
    double mu_prev = (mu_post * nu_post - mean * num_data) / nu_prev;
    double two_alpha_prev = two_alpha_post - num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_prev = beta_post - 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_prev, two_alpha_prev, beta_prev,
                                                          fctr->dp->hdp->log_sum_memo);

    cache_base_factor_params(fctr, mu_prev, nu_prev, two_alpha_prev, beta_prev, log_post_term);
}

double factor_parent_joint_log_likelihood(Factor* fctr, Factor* parent) {
    Factor* base_fctr = get_base_factor(parent);
    DirichletProcess* dp = fctr->dp;

    double num_reassign = (double) dp->cached_factor_size;
    double mean_reassign = dp->cached_factor_mean;
    double sum_sq_devs = dp->cached_factor_sum_sq_dev;

    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double nu_numer = nu_denom + num_reassign;

    double mean_dev = mean_reassign - mu_denom;
    double sq_mean_dev = nu_denom * num_reassign * mean_dev * mean_dev / nu_numer;

    double two_alpha_numer = two_alpha_denom + num_reassign;
    double beta_numer = beta_denom + 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_denom = param_array[4];
    double log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                      dp->hdp->log_sum_memo);

    return -0.5 * num_reassign * log(2.0 * M_PI) + log_numer - log_denom;
}

double data_pt_factor_parent_likelihood(Factor* data_pt_fctr, Factor* parent) {
    if (data_pt_fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Can only access data point likelihood for data point factors.\n");
        exit(EXIT_FAILURE);
    }

    double data_pt = get_factor_data_pt(data_pt_fctr);
    Factor* base_fctr = get_base_factor(parent);
    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double nu_numer = nu_denom + 1.0;

    double mean_dev = data_pt - mu_denom;
    double sq_mean_dev = nu_denom * mean_dev * mean_dev / nu_numer;

    double two_alpha_numer = two_alpha_denom + 1.0;
    double beta_numer = beta_denom + 0.5 * sq_mean_dev;

    double log_denom = param_array[4];
    double log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                      base_fctr->dp->hdp->log_sum_memo);

    return (1.0 / sqrt(2.0 * M_PI)) * exp(log_numer - log_denom);
}

void evaluate_posterior_predictive(Factor* base_fctr, double* x, double* pdf_out, int length,
							       SumOfLogsMemo* log_sum_memo) {
    if (base_fctr->factor_type != BASE) {
        fprintf(stderr, "Can only evaluate posterior predictive of base factors.\n");
        exit(EXIT_FAILURE);
    }

    double* param_array = base_fctr->factor_data;

    double mu_denom = param_array[0];
    double nu_denom = param_array[1];
    double two_alpha_denom = param_array[2];
    double beta_denom = param_array[3];

    double log_denom = param_array[4];

    double nu_numer = nu_denom + 1.0;
    double two_alpha_numer = two_alpha_denom + 1.0;
    double nu_ratio = nu_denom / nu_numer;
    double pi_factor = (1.0 / sqrt(2.0 * M_PI));

    double mean_dev;
    double sq_mean_dev;
    double beta_numer;
    double log_numer;
    for (int i = 0; i < length; i++) {
        mean_dev = x[i] - mu_denom;
        sq_mean_dev = nu_ratio * mean_dev * mean_dev;
        beta_numer = beta_denom + 0.5 * sq_mean_dev;

        log_numer = log_posterior_conditional_term(nu_numer, two_alpha_numer, beta_numer,
                                                   log_sum_memo);

        pdf_out[i] = pi_factor * exp(log_numer - log_denom);
    }
}

void evaluate_prior_predictive(HierarchicalDirichletProcess* hdp,
                               double* x, double* pdf_out, int length) {
    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    int two_alpha = (int) hdp->two_alpha;
    double beta = hdp->beta;

    double nu_factor = nu / (2.0 * (nu + 1.0) * beta);
    double alpha_term = exp(log_gamma_half(two_alpha + 1, hdp->log_sum_memo)
                            - log_gamma_half(two_alpha, hdp->log_sum_memo));
    double beta_term = sqrt(nu_factor / M_PI);
    double constant_term = alpha_term * beta_term;
    double power = -0.5 * (two_alpha + 1.0);

    for (int i = 0; i < length; i++) {
        double dev = x[i] - mu;
        double var_term = pow(1.0 + nu_factor * dev * dev, power);

        pdf_out[i] = constant_term * var_term;
    }
}

double prior_likelihood(HierarchicalDirichletProcess* hdp, Factor* fctr) {
    if (fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Cannot calculate point prior likelihood from non-data point factor.\n");
    }

    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    double dbl_two_alpha = hdp->two_alpha;
    int two_alpha = (int) dbl_two_alpha;
    double beta = hdp->beta;

    double data_pt = get_factor_data_pt(fctr);
    double dev = data_pt - mu;

    double alpha_term = exp(log_gamma_half(two_alpha + 1, hdp->log_sum_memo)
                        - log_gamma_half(two_alpha, hdp->log_sum_memo));
    double nu_term = nu / (2.0 * (nu + 1.0) * beta);
    double beta_term = pow(1.0 + nu_term * dev * dev, -0.5 * (dbl_two_alpha + 1.0));
    
    return alpha_term * sqrt(nu_term / M_PI) * beta_term;
}

double prior_joint_log_likelihood(HierarchicalDirichletProcess* hdp, Factor* fctr) {
    if (fctr->factor_type != MIDDLE) {
        fprintf(stderr, "Cannot calculate joint prior likelihood from non-middle factor.\n");
    }
    
    double mu = hdp->mu;
    double nu = hdp->nu;
    double dbl_two_alpha = hdp->two_alpha;
    int two_alpha = (int) dbl_two_alpha;
    double beta = hdp->beta;
    
    DirichletProcess* dp = fctr->dp;
    int num_reassign = dp->cached_factor_size;
    double dbl_reassign = (double) num_reassign;
    double mean_reassign = dp->cached_factor_mean;
    double sum_sq_devs = dp->cached_factor_sum_sq_dev;
    
    double mean_dev = mean_reassign - mu;
    double sq_mean_dev = nu * dbl_reassign * mean_dev * mean_dev / (nu + dbl_reassign);
    
    double log_alpha_term = log_gamma_half(two_alpha + num_reassign, hdp->log_sum_memo)
                             - log_gamma_half(two_alpha, hdp->log_sum_memo);
    double log_nu_term = 0.5 * (log(nu) - log(nu + dbl_reassign));
    double log_beta_term_1 = dbl_two_alpha * log(beta);
    double log_beta_term_2 = (dbl_two_alpha + dbl_reassign)
                              * log(beta + 0.5 * (sum_sq_devs + sq_mean_dev));
    return log_alpha_term + log_nu_term + 0.5 * (log_beta_term_1 - log_beta_term_2);
}

double unobserved_factor_likelihood(Factor* fctr, DirichletProcess* dp) {
    DirichletProcess* parent_dp = dp->parent;
    if (parent_dp == NULL) {
        return prior_likelihood(dp->hdp, fctr);
    }
    else {
        double parent_gamma = *(parent_dp->gamma);
        double likelihood = 0.0;
        
        stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
        Factor* parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
        double fctr_size;
        while (parent_fctr != NULL) {
            fctr_size = (double) stSet_size(parent_fctr->children);
            likelihood += fctr_size * data_pt_factor_parent_likelihood(fctr, parent_fctr);
            parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
        }
        stSet_destructIterator(parent_fctr_iter);
        
        likelihood += parent_gamma * unobserved_factor_likelihood(fctr, parent_dp);
        
        likelihood /= (parent_gamma + (double) parent_dp->num_factor_children);
        
        return likelihood;
    }
}

double unobserved_factor_joint_log_likelihood(Factor* fctr, DirichletProcess* dp) {
    DirichletProcess* parent_dp = dp->parent;
    if (parent_dp == NULL) {
        return prior_joint_log_likelihood(dp->hdp, fctr);
    }
    else {
        double parent_gamma = *(parent_dp->gamma);
        
        double log_likelihood = MINUS_INF;
        stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
        Factor* parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
        double log_fctr_size;
        while (parent_fctr != NULL) {
            log_fctr_size = log((double) stSet_size(parent_fctr->children));
            log_likelihood = add_logs(log_likelihood,
                                      log_fctr_size + factor_parent_joint_log_likelihood(fctr, parent_fctr));
            parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
        }
        stSet_destructIterator(parent_fctr_iter);
        
        log_likelihood = add_logs(log_likelihood,
                                  log(parent_gamma) + unobserved_factor_joint_log_likelihood(fctr, parent_dp));
        
        log_likelihood -= log(parent_gamma + (double) parent_dp->num_factor_children);
        
        return log_likelihood;
    }
}

DirichletProcess* new_dir_proc() {
    DirichletProcess* dp = (DirichletProcess*) malloc(sizeof(DirichletProcess));

    dp->gamma = NULL;
    dp->depth = 0;
    dp->parent = NULL;
    dp->children = stList_construct();
    dp->factors = stSet_construct();
    dp->num_factor_children = 0;

    dp->cached_factor_mean = 0.0;
    dp->cached_factor_sum_sq_dev = 0.0;
    dp->cached_factor_size = 0;

    dp->base_factor_wt = 0.0;
    dp->posterior_predictive = NULL;
    dp->spline_slopes = NULL;

    dp->observed = false;
    return dp;
}

void clear_factor_tree(Factor* fctr) {
    if (fctr->children != NULL) {
        stSetIterator* child_fctr_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            clear_factor_tree(child_fctr);
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
    }
    else {
        // note: this will trigger automatic destruction of parent factors
        destroy_factor(fctr);
    }
}

void destroy_dir_proc_factor_tree(DirichletProcess* dp) {
    if (stSet_size(dp->factors) == 0) {
        return;
    }
    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    while (fctr != NULL) {
        clear_factor_tree(fctr);
        fctr = (Factor*) stSet_getNext(fctr_iter);
    }
    stSet_destructIterator(fctr_iter);
}

void destroy_dir_proc(DirichletProcess* dp) {
    destroy_dir_proc_factor_tree(dp);
    stSet_destruct(dp->factors);
    
    if (dp->children != NULL) {
        stListIterator* st_iterator = stList_getIterator(dp->children);
        DirichletProcess* dp_child = (DirichletProcess*) stList_getNext(st_iterator);
        while (dp_child != NULL) {
            destroy_dir_proc(dp_child);
            dp_child = (DirichletProcess*) stList_getNext(st_iterator);
        }
        stList_destructIterator(st_iterator);
        
        stList_destruct(dp->children);
    }

    if (dp->parent != NULL) {
        stList_removeItem(dp->parent->children, (void*) dp);
    }

    free(dp->posterior_predictive);
    free(dp->spline_slopes);
    free(dp);
}

// fixed concentration parameters
HierarchicalDirichletProcess* new_hier_dir_proc(int num_dps, int depth, double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop, int sampling_grid_length, double mu,
                                                double nu, double alpha, double beta) {
    if (nu <= 0.0) {
        fprintf(stderr, "nu parameter of Normal-Inverse Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    //if (alpha <= 0.0) {
    //    fprintf(stderr, "alpha parameter of Normal-Inverse Gamma distribution must be positive.\n");
    //    exit(EXIT_FAILURE);
    //}
    if (beta <= 0.0) {
        fprintf(stderr, "beta parameter of Normal-Inverse Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (2 * alpha != (int) 2 * alpha || alpha <= 1.0) {
            fprintf(stderr, "Normal-Inverse Gamma parameter 'alpha' must be integer or half-integer > 1.0.\n");
            exit(EXIT_FAILURE);
        }
    if (gamma != NULL) {
        for (int i = 0; i < depth; i++) {
            if (gamma[i] <= 0) {
                fprintf(stderr, "Concentration parameter gamma must be postive.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    double* grid = linspace(sampling_grid_start, sampling_grid_stop, sampling_grid_length);

    HierarchicalDirichletProcess* hdp = (HierarchicalDirichletProcess*) malloc(sizeof(HierarchicalDirichletProcess));

    // normal-inverse gamma parameters

    hdp->mu = mu;
    hdp->nu = nu;
    hdp->two_alpha = 2.0 * alpha;
    hdp->beta = beta;

    hdp->gamma = gamma;
    hdp->depth = depth;

    hdp->finalized = false;

    hdp->num_dps = num_dps;
    DirichletProcess** dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    for (int i = 0; i < num_dps; i++) {
        DirichletProcess* dp = new_dir_proc();
        dp->id = i;
        dp->hdp = hdp;
        dps[i] = dp;
    }

    hdp->dps = dps;
    hdp->base_dp = NULL;

    hdp->sampling_grid = grid;
    hdp->grid_length = sampling_grid_length;
    hdp->samples_taken = 0;
    hdp->splines_finalized = false;

    hdp->data = NULL;
    hdp->data_pt_dp_id = NULL;
    hdp->data_length = 0;

    hdp->log_sum_memo = new_log_sum_memo();

    hdp->sample_gamma = false;
    hdp->s_aux_vector = NULL;
    hdp->w_aux_vector = NULL;

    return hdp;
}

// Gamma prior on concentration parameters
HierarchicalDirichletProcess* new_hier_dir_proc_2(int num_dps, int depth, double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start, double sampling_grid_stop,
                                                  int sampling_grid_length, double mu, double nu, double alpha,
                                                  double beta) {

    for (int i = 0; i < depth; i++) {
        if (gamma_alpha[i] <= 0.0) {
            fprintf(stderr, "alpha parameter of Gamma distribution must be positive.\n");
            exit(EXIT_FAILURE);
        }
        if (gamma_beta[i] <= 0.0) {
            fprintf(stderr, "beta parameter of Gamma distribution must be positive.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dps, depth, NULL, sampling_grid_start, sampling_grid_stop,
                                                          sampling_grid_length, mu, nu, alpha, beta);

    hdp->sample_gamma = true;
    
    hdp->gamma_alpha = gamma_alpha;
    hdp->gamma_beta = gamma_beta;

    double* w = (double*) malloc(sizeof(double) * num_dps);
    hdp->w_aux_vector = w;
    bool* s = (bool*) malloc(sizeof(bool) * num_dps);
    hdp->s_aux_vector = s;

    for (int i = 0; i < num_dps; i++) {
        w[i] = 1.0;
        s[i] = false;
    }
    
    // init to prior expected value
    double* gamma = (double*) malloc(sizeof(double) * depth);
    hdp->gamma = gamma;
    for (int i = 0; i < depth; i++) {
        gamma[i] = gamma_alpha[i] / gamma_beta[i];
    }
    
    
    return hdp;
}

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp) {
    destroy_dir_proc(hdp->base_dp);
    free(hdp->gamma);
    free(hdp->data);
    free(hdp->data_pt_dp_id);
    free(hdp->dps);
    free(hdp->sampling_grid);
    destroy_log_sum_memo(hdp->log_sum_memo);
    free(hdp->gamma_alpha);
    free(hdp->gamma_beta);
    free(hdp->w_aux_vector);
    free(hdp->s_aux_vector);
    free(hdp);
}

void establish_base_dp(HierarchicalDirichletProcess* hdp) {
    DirichletProcess** dps = hdp->dps;
    int num_dps = hdp->num_dps;

    DirichletProcess* dp;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];
        if (dp->parent == NULL) {
            if (hdp->base_dp == NULL) {
                hdp->base_dp = dp;
            }
            else {
                fprintf(stderr, "Hierarchical Dirichlet process contains orphaned Dirichlet process.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (hdp->base_dp == NULL) {
        fprintf(stderr, "Hierarchical Dirichlet process does not contain base Dirichlet process.\n");
        exit(EXIT_FAILURE);
    }
}

// DFS to verify that Dirichlet processes follow tree structure
void verify_dp_tree(HierarchicalDirichletProcess* hdp) {
    int num_dps = hdp->num_dps;
    bool* visited = (bool*) malloc(sizeof(bool) * num_dps);
    for (int i = 0; i < num_dps; i++) {
        visited[i] = false;
    }

    DirichletProcess* base_dp = hdp->base_dp;
    stList* stck = stList_construct();
    stList_append(stck, (void*) base_dp);

    DirichletProcess* dp;
    while (stList_length(stck) > 0) {
        dp = (DirichletProcess*) stList_pop(stck);
        if (visited[dp->id]) {
            fprintf(stderr, "Hierarchical Dirichlet process does not have tree structure.\n");
            exit(EXIT_FAILURE);
        }
        visited[dp->id] = true;

        stListIterator* child_iter = stList_getIterator(dp->children);
        DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
        while (child != NULL) {
            stList_append(stck, (void*) child);
            child = (DirichletProcess*) stList_getNext(child_iter);
        }
        stList_destructIterator(child_iter);
    }
    stList_destruct(stck);
    free(visited);
}

void verify_tree_depth(HierarchicalDirichletProcess* hdp, DirichletProcess* dp, int current_depth,
                       int leaf_depth) {
    dp->gamma = &(hdp->gamma[current_depth]);
    dp->depth = current_depth;

    if (stList_length(dp->children) == 0) {
        if (current_depth != leaf_depth) {
            fprintf(stderr, "Hierarchical Dirichlet process has leaf Dirichlet process at incorrect depth.\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        stListIterator* st_iterator = stList_getIterator(dp->children);
        DirichletProcess* child = (DirichletProcess*) stList_getNext(st_iterator);
        while (child != NULL) {
            verify_tree_depth(hdp, child, current_depth + 1, leaf_depth);
            child = (DirichletProcess*) stList_getNext(st_iterator);
        }
        stList_destructIterator(st_iterator);
    }
}

void verify_valid_dp_assignments(HierarchicalDirichletProcess* hdp) {

    int length = hdp->data_length;
    int num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    int* dp_ids = hdp->data_pt_dp_id;

    int id;
    DirichletProcess* dp;
    for (int i = 0; i < length; i++) {
        id = dp_ids[i];
        if (id >= num_dps || id < 0) {
            fprintf(stderr, "Data point is assigned to non-existent Dirichlet process.\n");
            exit(EXIT_FAILURE);
        }

        dp = dps[id];
        if (stList_length(dp->children) > 0) {
            fprintf(stderr, "Data point cannot be assigned to non-leaf Dirichlet process.\n");
            exit(EXIT_FAILURE);
        }
    }
}

void mark_observed_dps(HierarchicalDirichletProcess* hdp) {
    // mark newly observed dps
    int length = hdp->data_length;
    DirichletProcess** dps = hdp->dps;
    int* dp_ids = hdp->data_pt_dp_id;
    int grid_length = hdp->grid_length;

    DirichletProcess* dp;
    double* pdf;
    int id;
    for (int i = 0; i < length; i++) {
        id = dp_ids[i];
        dp = dps[id];
        while (dp != NULL) {
            if (dp->observed) {
                break;
            }
            dp->observed = true;

            pdf = (double*) malloc(sizeof(double) * grid_length);
            dp->posterior_predictive = pdf;
            for (int j = 0; j < grid_length; j++) {
                pdf[j] = 0.0;
            }

            dp = dp->parent;
        }
    }
}

void init_factors_internal(DirichletProcess* dp, Factor* parent_fctr, stList** data_pt_fctr_lists) {
    if (!dp->observed) {
        return;
    }
    Factor* fctr = new_middle_factor(dp);
    fctr->parent = parent_fctr;
    stSet_insert(parent_fctr->children, (void*) fctr);
    
    if (stList_length(dp->children) == 0) {
        stSet* children = fctr->children;
        stListIterator* data_pt_fctr_iter = stList_getIterator(data_pt_fctr_lists[dp->id]);
        Factor* data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        while (data_pt_fctr != NULL) {
            data_pt_fctr->parent = fctr;
            stSet_insert(children, (void*) data_pt_fctr);
            data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        }
        stList_destructIterator(data_pt_fctr_iter);
    }
    else {
        stListIterator* child_dp_iter = stList_getIterator(dp->children);
        DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        while (child_dp != NULL) {
            init_factors_internal(child_dp, fctr, data_pt_fctr_lists);
            child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        }
        stList_destructIterator(child_dp_iter);
    }
}

void init_factors(HierarchicalDirichletProcess* hdp) {
    int data_length = hdp->data_length;
    int* data_pt_dp_id = hdp->data_pt_dp_id;
    DirichletProcess** dps = hdp->dps;
    int num_dps = hdp->num_dps;
    
    stList** data_pt_fctr_lists = (stList**) malloc(sizeof(stList*) * num_dps);
    for (int i = 0; i < num_dps; i++) {
        data_pt_fctr_lists[i] = NULL;
    }
    
    Factor* data_pt_fctr;
    int dp_id;
    stList* fctr_list;
    for (int data_pt_idx = 0; data_pt_idx < data_length; data_pt_idx++) {
        dp_id = data_pt_dp_id[data_pt_idx];
        fctr_list = data_pt_fctr_lists[dp_id];
        if (fctr_list == NULL) {
            fctr_list = stList_construct();
            data_pt_fctr_lists[dp_id] = fctr_list;
        }
        data_pt_fctr = new_data_pt_factor(hdp, data_pt_idx);
        stList_append(fctr_list, (void*) data_pt_fctr);
    }
    
    DirichletProcess* base_dp = hdp->base_dp;
    Factor* root_factor = new_base_factor(hdp);
    
    stListIterator* child_dp_iter = stList_getIterator(base_dp->children);
    DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    while (child_dp != NULL) {
        init_factors_internal(child_dp, root_factor, data_pt_fctr_lists);
        child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    }
    stList_destructIterator(child_dp_iter);
    
    for (int i = 0; i < num_dps; i++) {
        if (data_pt_fctr_lists[i] != NULL) {
            stList_destruct(data_pt_fctr_lists[i]);
        }
    }
    free(data_pt_fctr_lists);
    
    double mean, sum_sq_devs;
    int num_data;
    get_factor_stats(root_factor, &mean, &sum_sq_devs, &num_data);
    add_update_base_factor_params(root_factor, mean, sum_sq_devs, (double) num_data);
    
    int fctr_child_count;
    DirichletProcess* dp;
    Factor* fctr;
    stSetIterator* fctr_iter;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];
        
        fctr_child_count = 0;

        fctr_iter = stSet_getIterator(dp->factors);
        fctr = (Factor*) stSet_getNext(fctr_iter);
        while (fctr != NULL) {
            fctr_child_count += stSet_size(fctr->children);
            fctr = (Factor*) stSet_getNext(fctr_iter);
        }
        stSet_destructIterator(fctr_iter);

        dp->num_factor_children = fctr_child_count;
    }
}

void finalize_data(HierarchicalDirichletProcess* hdp) {
    verify_valid_dp_assignments(hdp);
    mark_observed_dps(hdp);
    init_factors(hdp);
}

void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int child_id, int parent_id) {
    if (hdp->finalized) {
        fprintf(stderr, "Hierarchical Dirichlet process structure has been finalized. Cannot set new parent.\n");
        exit(EXIT_FAILURE);
    }

    if (child_id >= hdp->num_dps || parent_id >= hdp->num_dps || child_id < 0 || parent_id < 0) {
        fprintf(stderr, "Dirichlet process ID does not exist.\n");
        exit(EXIT_FAILURE);
    }

    DirichletProcess* child_dp = hdp->dps[child_id];
    DirichletProcess* parent_dp = hdp->dps[parent_id];

    if (child_dp->parent != NULL) {
        fprintf(stderr, "Dirichlet process already has parent.\n");
        exit(EXIT_FAILURE);
    }

    child_dp->parent = parent_dp;
    stList_append(parent_dp->children, (void*) child_dp);
}

void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int* dp_ids, int length) {
    if (hdp->data != NULL) {
        fprintf(stderr, "Hierarchical Dirichlet process must be reset before passing new data.\n");
        exit(EXIT_FAILURE);
    }

    hdp->data = data;
    hdp->data_pt_dp_id = dp_ids;
    hdp->data_length = length;

    if (hdp->finalized) {
        finalize_data(hdp);
    }
}

void finalize_hdp_structure(HierarchicalDirichletProcess* hdp) {
    establish_base_dp(hdp);
    verify_dp_tree(hdp);
    verify_tree_depth(hdp, hdp->base_dp, 0, hdp->depth - 1);

    if (hdp->data != NULL) {
        finalize_data(hdp);
    }

    hdp->finalized = true;
}

void reset_hdp_data(HierarchicalDirichletProcess* hdp) {
    if (hdp->data == NULL && hdp->data_pt_dp_id == NULL) {
        return;
    }
    
    free(hdp->data);
    hdp->data = NULL;

    free(hdp->data_pt_dp_id);
    hdp->data_pt_dp_id = NULL;

    DirichletProcess** dps = hdp->dps;
    int num_dps = hdp->num_dps;
    
    destroy_dir_proc_factor_tree(hdp->base_dp);
    
    DirichletProcess* dp;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];

        free(dp->posterior_predictive);
        dp->posterior_predictive = NULL;

        free(dp->spline_slopes);
        dp->spline_slopes = NULL;

        dp->observed = false;
    }
    
    hdp->splines_finalized = false;
    
    hdp->samples_taken = 0;
    
    if (hdp->sample_gamma) {
        double* gamma = hdp->gamma;
        double* gamma_alpha = hdp->gamma_alpha;
        double* gamma_beta = hdp->gamma_beta;
        
        for (int depth = 0; depth < hdp->depth; depth++) {
            gamma[depth] = gamma_alpha[depth] / gamma_beta[depth];
        }
        
        double* w = hdp->w_aux_vector;
        bool* s = hdp->s_aux_vector;
        
        for (int i = 0; i < num_dps; i++) {
            w[i] = 1.0;
            s[i] = false;
        }
    }
}

void unassign_from_parent(Factor* fctr) {
    if (fctr->factor_type == BASE) {
        fprintf(stderr, "Cannot unassign base factor's parent.\n");
        exit(EXIT_FAILURE);
    }
    Factor* parent = fctr->parent;
    Factor* base_fctr = get_base_factor(parent);
    DirichletProcess* base_dp = base_fctr->dp;
    
    stSet_remove(parent->children, (void*) fctr);
    fctr->parent = NULL;
    (parent->dp->num_factor_children)--;
    
    if (stSet_size(parent->children) == 0) {
        destroy_factor(parent);
    }

    int num_reassign;
    double mean_reassign;
    double sum_sq_devs;
    
    get_factor_stats(fctr, &mean_reassign, &sum_sq_devs, &num_reassign);
    
    // check to see if base factor has been destroyed
    if (stSet_search(base_dp->factors, (void*) base_fctr) != NULL) {
        remove_update_base_factor_params(base_fctr, mean_reassign, sum_sq_devs, (double) num_reassign);
    }
    
    DirichletProcess* dp = fctr->dp;
    if (dp != NULL) {
        dp->cached_factor_mean = mean_reassign;
        dp->cached_factor_size = num_reassign;
        dp->cached_factor_sum_sq_dev = sum_sq_devs;
    }
}

void assign_to_parent(Factor* fctr, Factor* parent) {
    if (fctr->factor_type == BASE) {
        fprintf(stderr, "Cannot assign base factor to a parent.\n");
        exit(EXIT_FAILURE);
    }
    
    if (parent->factor_type == DATA_PT) {
        fprintf(stderr, "Cannot assign data point factor to be parent.\n");
        exit(EXIT_FAILURE);
    }
    
    fctr->parent = parent;
    stSet_insert(parent->children, (void*) fctr);
    (parent->dp->num_factor_children)++;

    Factor* base_fctr = get_base_factor(parent);
    if (base_fctr == NULL) {
        return;
    }

    if (fctr->factor_type == DATA_PT) {
        double data_pt = get_factor_data_pt(fctr);
        add_update_base_factor_params(base_fctr, data_pt, 0.0, 1.0);
    }
    else {
        DirichletProcess* dp = fctr->dp;
        add_update_base_factor_params(base_fctr, dp->cached_factor_mean, dp->cached_factor_sum_sq_dev,
                                      (double) dp->cached_factor_size);
    }
}

Factor* sample_from_data_pt_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != DATA_PT) {
        fprintf(stderr, "Attempted a data point factor sample from non-data point factor.\n");
        exit(EXIT_FAILURE);
    }
    
    stSet* pool = dp->factors;
    int num_fctrs = stSet_size(pool);
    
    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);

    double* cdf = (double*) malloc(sizeof(double) * (num_fctrs + 1));
    double cumul = 0.0;

    stSetIterator* pool_iter = stSet_getIterator(pool);
    Factor* fctr_option;
    double fctr_size;
    for (int i = 0; i < num_fctrs; i++) {
        fctr_option = (Factor*) stSet_getNext(pool_iter);
        fctr_order[i] = fctr_option;
        
        fctr_size = (double) stSet_size(fctr_option->children);
        cumul += fctr_size * data_pt_factor_parent_likelihood(fctr, fctr_option);
        cdf[i] = cumul;
    }
    stSet_destructIterator(pool_iter);
    
    double gamma_param = *(dp->gamma);
    cumul += gamma_param * unobserved_factor_likelihood(fctr, dp);
    cdf[num_fctrs] = cumul;
    
    int choice_idx = bisect_left(rand_uniform(cumul), cdf, num_fctrs + 1);
    free(cdf);
    
    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_data_pt_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
    }

    free(fctr_order);

    return fctr_choice;
}

Factor* sample_from_middle_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != MIDDLE) {
        fprintf(stderr, "Attempted a middle factor sample from non-middle factor.\n");
        exit(EXIT_FAILURE);
    }
    
    stSet* pool = dp->factors;
    int num_fctrs = stSet_size(pool);
    int num_choices = num_fctrs + 1;

    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);
    double* log_probs = (double*) malloc(sizeof(double) * num_choices);
    
    stSetIterator* pool_iter = stSet_getIterator(pool);
    Factor* fctr_option;
    for (int i = 0; i < num_fctrs; i++) {
        fctr_option = (Factor*) stSet_getNext(pool_iter);
        fctr_order[i] = fctr_option;
        log_probs[i] = factor_parent_joint_log_likelihood(fctr, fctr_option);
    }
    stSet_destructIterator(pool_iter);
    
    log_probs[num_fctrs] = unobserved_factor_joint_log_likelihood(fctr, dp);
    
    double* cdf = (double*) malloc(sizeof(double) * num_choices);
    double cumul = 0.0;
    double normalizing_const = max(log_probs, num_choices);
    
    double fctr_size;
    for (int i = 0; i < num_fctrs; i++) {
        fctr_size = (double) stSet_size(fctr_option->children);
        cumul += fctr_size * exp(log_probs[i] - normalizing_const);
        cdf[i] = cumul;
    }
    
    double gamma_param = *(dp->gamma);
    cumul += gamma_param * exp(log_probs[num_fctrs] - normalizing_const);
    cdf[num_fctrs] = cumul;
    
    free(log_probs);

    int choice_idx = bisect_left(rand_uniform(cumul), cdf, num_choices);
    free(cdf);

    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_middle_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
    }
    free(fctr_order);

    return fctr_choice;
}

Factor* sample_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type == DATA_PT) {
        return sample_from_data_pt_factor(fctr, dp);
    }
    else if (fctr->factor_type == MIDDLE) {
        return sample_from_middle_factor(fctr, dp);
    }
    else {
        fprintf(stderr, "Cannot sample base factor parent assignments.\n");
        exit(EXIT_FAILURE);
    }
}

void gibbs_factor_iteration(Factor* fctr) {
    DirichletProcess* parent_dp = fctr->parent->dp;
    unassign_from_parent(fctr);
    Factor* new_parent = sample_factor(fctr, parent_dp);
    assign_to_parent(fctr, new_parent);
}

void cache_prior_contribution(DirichletProcess* dp, double parent_prior_prod) {
    if (!(dp->observed)) {
        return;
    }
    double gamma_param = *(dp->gamma);
    double total_children = (double) dp->num_factor_children;
    double prior_prod = (gamma_param / (gamma_param + total_children)) * parent_prior_prod;
    dp->base_factor_wt += prior_prod;

    stListIterator* child_iter = stList_getIterator(dp->children);
    DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
    while (child != NULL) {
        cache_prior_contribution(child, prior_prod);
        child = (DirichletProcess*) stList_getNext(child_iter);
    }
    stList_destructIterator(child_iter);
}

void cache_base_factor_weight(Factor* fctr) {
    DirichletProcess* dp = fctr->dp;
    
    double gamma_param = *(dp->gamma);
    double total_children = (double) dp->num_factor_children;
    double wt = ((double) stSet_size(fctr->children)) / (gamma_param + total_children);
    dp->base_factor_wt += wt;
    
    if (stList_length(dp->children) > 0) {
        stSetIterator* child_fctr_iter = stSet_getIterator(fctr->children);
        Factor* child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            cache_base_factor_weight(child_fctr);
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        
        stListIterator* child_dp_iter = stList_getIterator(dp->children);
        DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        while (child_dp != NULL) {
            cache_prior_contribution(child_dp, wt);
            child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        }
        stList_destructIterator(child_dp_iter);
    }
}

void push_factor_distr(DirichletProcess* dp, double* distr, int length) {
    double* sample_collector = dp->posterior_predictive;
    double wt = dp->base_factor_wt;

    for (int i = 0; i < length; i++) {
        sample_collector[i] += wt * distr[i];
    }

    dp->base_factor_wt = 0.0;

    stListIterator* child_iter = stList_getIterator(dp->children);
    DirichletProcess* child = (DirichletProcess*) stList_getNext(child_iter);
    while (child != NULL) {
        if (child->observed) {
            push_factor_distr(child, distr, length);
        }
        child = (DirichletProcess*) stList_getNext(child_iter);
    }
    stList_destructIterator(child_iter);
}

void take_distr_sample(HierarchicalDirichletProcess* hdp) {
    DirichletProcess* base_dp = hdp->base_dp;

    double* grid = hdp->sampling_grid;
    int length = hdp->grid_length;
    double* pdf = (double*) malloc(sizeof(double) * length);
    
    SumOfLogsMemo* log_sum_memo = hdp->log_sum_memo;

    stSetIterator* base_fctr_iter = stSet_getIterator(base_dp->factors);
    Factor* base_fctr = (Factor*) stSet_getNext(base_fctr_iter);
    while (base_fctr != NULL) {
        cache_base_factor_weight(base_fctr);
        evaluate_posterior_predictive(base_fctr, grid, pdf, length, log_sum_memo);
        push_factor_distr(base_dp, pdf, length);

        base_fctr = (Factor*) stSet_getNext(base_fctr_iter);
    }
    stSet_destructIterator(base_fctr_iter);
    
    cache_prior_contribution(base_dp, 1.0);
    evaluate_prior_predictive(hdp, grid, pdf, length);
    push_factor_distr(base_dp, pdf, length);

    (hdp->samples_taken)++;
    
    free(pdf);
}

// Knuth shuffle algorithm
DirichletProcess** get_shuffled_dps(HierarchicalDirichletProcess* hdp) {
    int num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess** shuffled_dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    int pos;
    for (int i = 0; i < num_dps; i++) {
        pos = rand() % (i + 1);
        shuffled_dps[i] = shuffled_dps[pos];
        shuffled_dps[pos] = dps[i];
    }
    return shuffled_dps;
}

void sample_dp_factors(DirichletProcess* dp, int* iter_counter, int burn_in, int thinning,
                       int* sample_counter, int num_samples) {
    
    if (!dp->observed) {
        return;
    }

    int iter = *iter_counter;
    int samples_taken = *sample_counter;
    
    // have to pre-allocate the array of sampling factors in case reassignment triggers
    // destruction of the set the iterator is iterating through
    int num_factor_children = dp->num_factor_children;
    Factor** sampling_fctrs = (Factor**) malloc(sizeof(Factor*) * num_factor_children);
    int i = 0;
    
    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    stSetIterator* child_fctr_iter;
    Factor* child_fctr;
    while (fctr != NULL) {
        child_fctr_iter = stSet_getIterator(fctr->children);
        child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            sampling_fctrs[i] = child_fctr;
            i++;
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        fctr = (Factor*) stSet_getNext(fctr_iter);
    }
    stSet_destructIterator(fctr_iter);
    
    for (int j = 0; j < num_factor_children; j++) {
        gibbs_factor_iteration(sampling_fctrs[j]);
        iter++;
        
        if (iter % thinning == 0) {
            
            if (iter > burn_in) {
                take_distr_sample(dp->hdp);
                samples_taken++;
                
                if (samples_taken >= num_samples) {
                    break;
                }
            }
        }
    }
    free(sampling_fctrs);
    
    *sample_counter = samples_taken;
    *iter_counter = iter;
}

double sample_auxilliary_w(DirichletProcess* dp) {
    return (double) genbet((float) *(dp->gamma) + 1.0, (float) dp->num_factor_children);
}

bool sample_auxilliary_s(DirichletProcess* dp) {
    double num_children = (double) dp->num_factor_children;
    return rand_bernoulli(num_children / (num_children + *(dp->gamma)));
}

void sample_gamma_aux_vars(HierarchicalDirichletProcess* hdp) {
    double* w = hdp->w_aux_vector;
    bool* s = hdp->s_aux_vector;

    DirichletProcess** dps = hdp->dps;
    int num_dps = hdp->num_dps;
    DirichletProcess* dp;
    for (int id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        w[id] = sample_auxilliary_w(dp);
        s[id] = sample_auxilliary_s(dp);
    }
}

void sample_base_gamma_internal(HierarchicalDirichletProcess* hdp, double log_w, int num_factors) {
    // Escobar and West's (1995) algorithm
    DirichletProcess* base_dp = hdp->base_dp;
    double gamma_alpha = hdp->gamma_alpha[0];
    double gamma_beta = hdp->gamma_beta[0];

    double num_children = (double) base_dp->num_factor_children;

    double gamma_beta_post = gamma_beta - log_w;
    double gamma_alpha_post = gamma_alpha + (double) num_factors;

    double frac = (gamma_alpha_post - 1.0)
                  / (num_children * gamma_beta_post);

    double wt = frac / (1.0 + frac);
    // note: different parameterization switches alpha and beta
    float sample_gamma = wt * gengam(gamma_beta_post, gamma_alpha_post)
                         + (1 - wt) * gengam(gamma_beta_post, gamma_alpha_post - 1.0);

    hdp->gamma[0] = (double) sample_gamma;
}

void sample_middle_gammas_internal(HierarchicalDirichletProcess* hdp, int depth,
                         double sum_log_w, int sum_s, int num_depth_fctrs) {
    double gamma_alpha = hdp->gamma_alpha[depth];
    double gamma_beta = hdp->gamma_beta[depth];

    float gamma_alpha_post = (float) (gamma_alpha + (double) (num_depth_fctrs - sum_s));
    float gamma_beta_post = (float) (gamma_beta - sum_log_w);
    // note: different parameterization switches alpha and beta
    hdp->gamma[depth] = (double) gengam(gamma_beta_post, gamma_alpha_post);
}

void sample_gammas(HierarchicalDirichletProcess* hdp, int* iter_counter, int burn_in,
                   int thinning, int* sample_counter, int num_samples) {
    int iter = *iter_counter;
    int samples_taken = *sample_counter;

    int tree_depth = hdp->depth;
    double* w = hdp->w_aux_vector;
    bool* s = hdp->s_aux_vector;

    int* num_depth_fctrs = (int*) malloc(sizeof(int) * tree_depth);
    double* sum_log_w = (double*) malloc(sizeof(double) * tree_depth);
    int* sum_s = (int*) malloc(sizeof(int) * tree_depth);

    for (int depth = 0; depth < tree_depth; depth++) {
        num_depth_fctrs[depth] = 0;
        sum_log_w[depth] = 0.0;
        sum_s[depth] = 0;
    }

    int num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    int dp_depth;
    for (int id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        dp_depth = dp->depth;
        num_depth_fctrs[dp_depth] += stSet_size(dp->factors);
        sum_log_w[dp_depth] += log(w[id]);
        if (s[id]) sum_s[dp_depth]++;
    }

    for (int depth = 0; depth < tree_depth; depth++) {
        if (depth == 0) {
            sample_base_gamma_internal(hdp, sum_log_w[depth], num_depth_fctrs[depth]);
        }
        else {
            sample_middle_gammas_internal(hdp, depth, sum_log_w[depth],
                                          sum_s[depth], num_depth_fctrs[depth]);
        }
        iter++;

        if (iter % thinning == 0) {
            if (iter > burn_in) {
                take_distr_sample(dp->hdp);
                samples_taken++;
                
                if (samples_taken >= num_samples) {
                    break;
                }
            }
        }
    }
    free(sum_log_w);
    free(sum_s);
    free(num_depth_fctrs);

    *iter_counter = iter;
    *sample_counter = samples_taken;
}

void sample_gamma_params(HierarchicalDirichletProcess* hdp, int* iter_counter, int burn_in,
                         int thinning, int* sample_counter, int num_samples) {
    sample_gamma_aux_vars(hdp);
    sample_gammas(hdp, iter_counter, burn_in, thinning, sample_counter, num_samples);
}



int* snapshot_num_factors(HierarchicalDirichletProcess* hdp, int* length_out) {
    int length = hdp->num_dps;
    *length_out = length;
    int* snapshot = (int*) malloc(sizeof(int) * length);
    
    DirichletProcess** dps = hdp->dps;
    for (int i = 0; i < length; i++) {
        snapshot[i] = (int) stSet_size((dps[i])->factors);
    }
    return snapshot;
}

double* snapshot_gamma_params(HierarchicalDirichletProcess* hdp, int* length_out) {
    int length = hdp->depth;
    *length_out = length;
    
    double* snapshot = (double*) malloc(sizeof(double) * length);
    double* gammas = hdp->gamma;
    for (int i = 0; i < length; i++) {
        snapshot[i] = gammas[i];
    }
    return snapshot;
}

double snapshot_factor_log_likelihood(Factor* fctr) {
    
    double parent_prob;
    double cumul = 0.0;
    
    if (fctr->factor_type == BASE) {
        fprintf(stderr, "Cannot snapshot base factor log likelihood.\n");
        exit(EXIT_FAILURE);
    }
    else if (fctr->factor_type == DATA_PT) {
        Factor* parent_fctr = fctr->parent;
        DirichletProcess* parent_dp = parent_fctr->dp;
        stSet* pool = parent_dp->factors;
        int num_fctrs = stSet_size(pool);
        
        stSetIterator* pool_iter = stSet_getIterator(pool);
        Factor* fctr_option;
        double fctr_size;
        double prob;
        for (int i = 0; i < num_fctrs; i++) {
            fctr_option = (Factor*) stSet_getNext(pool_iter);
            fctr_size = (double) stSet_size(fctr_option->children);
            prob = fctr_size * data_pt_factor_parent_likelihood(fctr, fctr_option);
            cumul += prob;
            if (fctr_option == parent_fctr) {
                parent_prob = prob;
            }
        }
        stSet_destructIterator(pool_iter);
        
        double gamma_param = *(parent_dp->gamma);
        cumul += gamma_param * unobserved_factor_likelihood(fctr, parent_dp);
    }
    else {
        DirichletProcess* dp = fctr->dp;
        DirichletProcess* parent_dp = dp->parent;
        Factor* parent_fctr = fctr->parent;
        
        double mean, sum_sq_devs;
        int num_data;
        get_factor_stats(fctr, &mean, &sum_sq_devs, &num_data);
        
        dp->cached_factor_mean = mean;
        dp->cached_factor_size = num_data;
        dp->cached_factor_sum_sq_dev = sum_sq_devs;
        
        stSet* pool = parent_dp->factors;
        int num_fctrs = stSet_size(pool);
        int num_choices = num_fctrs + 1;
        
        double* log_probs = (double*) malloc(sizeof(double) * num_choices);
        
        stSetIterator* pool_iter = stSet_getIterator(pool);
        Factor* fctr_option;
        double log_prob;
        double parent_log_prob;
        double fctr_size;
        for (int i = 0; i < num_fctrs; i++) {
            fctr_option = (Factor*) stSet_getNext(pool_iter);
            fctr_size = (double) stSet_size(fctr_option->children);
            log_prob = factor_parent_joint_log_likelihood(fctr, fctr_option) + log(fctr_size);
            log_probs[i] = log_prob;
            if (fctr_option == parent_fctr) {
                parent_log_prob = log_prob;
            }
            
        }
        stSet_destructIterator(pool_iter);
        
        double gamma_param = *(dp->gamma);
        log_probs[num_fctrs] = unobserved_factor_joint_log_likelihood(fctr, parent_dp) + log(gamma_param);
        
        double normalizing_const = max(log_probs, num_choices);
        
        parent_prob = exp(parent_log_prob - normalizing_const);
        
        for (int i = 0; i < num_choices; i++) {
            cumul += exp(log_probs[i] - normalizing_const);;
        }
        
        free(log_probs);
    }
    
    // TODO: this is a hack, makes it inaccurate for the early iterations
    if (parent_prob == 0.0) {
        return 0.0;
    }
    return (log(parent_prob) - log(cumul)) / 1000.0;
}

double snapshot_dir_proc_log_likelihood(DirichletProcess* dp) {
    double log_likelihood = 0.0;
    
    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    
    stSetIterator* child_fctr_iter;
    Factor* child_fctr;
    while (fctr != NULL) {
        child_fctr_iter = stSet_getIterator(fctr->children);
        child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            log_likelihood += snapshot_factor_log_likelihood(child_fctr);
            child_fctr = (Factor*) stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        fctr = (Factor*) stSet_getNext(fctr_iter);
    }
    stSet_destructIterator(fctr_iter);
    return log_likelihood;
}

double snapshot_log_likelihood(HierarchicalDirichletProcess* hdp) {
    
    double log_likelihood = 0.0;
    
    int num_dps = hdp->num_dps;
    
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    for (int id = 0; id < num_dps; id++){
        dp = dps[id];
        
        if (!dp->observed) {
            continue;
        }
        
        log_likelihood += snapshot_dir_proc_log_likelihood(dp);
    }
    
    return log_likelihood;
}

void take_snapshot(HierarchicalDirichletProcess* hdp, int** num_dp_fctrs_out, int* num_dps_out,
                   double** gamma_params_out, int* num_gamma_params_out, double* log_likelihood_out) {
    
    *num_dp_fctrs_out = snapshot_num_factors(hdp, num_dps_out);
    *gamma_params_out = snapshot_gamma_params(hdp, num_gamma_params_out);
    *log_likelihood_out = snapshot_log_likelihood(hdp);
}

void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int num_samples, int burn_in,
                            int thinning) {
    
    execute_gibbs_sampling_with_snapshots(hdp, num_samples, burn_in, thinning, NULL, NULL);
}

void execute_gibbs_sampling_with_snapshots(HierarchicalDirichletProcess* hdp, int num_samples, int burn_in, int thinning,
                                           void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                           void* snapshot_func_args) {
    if (hdp->data == NULL || hdp->data_pt_dp_id == NULL) {
        fprintf(stderr, "Cannot perform Gibbs sampling before passing data to HDP.\n");
        exit(EXIT_FAILURE);
    }
    
    int iter_counter = 0;
    int sample_counter = 0;
    int num_dps = hdp->num_dps;

    DirichletProcess** sampling_dps;
    while (sample_counter < num_samples) {
        
        if (snapshot_func != NULL) {
            snapshot_func(hdp, snapshot_func_args);
        }
        
        sampling_dps = get_shuffled_dps(hdp);

        for (int i = 0; i < num_dps; i++) {
            sample_dp_factors(sampling_dps[i], &iter_counter, burn_in, thinning,
                              &sample_counter, num_samples);
            if (sample_counter >= num_samples) {
                break;
            }
        }

        free(sampling_dps);

        if (hdp->sample_gamma && sample_counter < num_samples) {
            sample_gamma_params(hdp, &iter_counter, burn_in, thinning, &sample_counter,
                                num_samples);
        }
    }
}

void finalize_distributions(HierarchicalDirichletProcess* hdp) {
    if (hdp->samples_taken <= 0) {
        fprintf(stderr, "Must perform Gibbs sampling before finalizing sampled distributions.\n");
        exit(EXIT_FAILURE);
    }
    
    if (hdp->splines_finalized) {
        fprintf(stderr, "Distributions have already been finalized.\n");
        exit(EXIT_FAILURE);
    }

    double inv_sample_size = 1.0 / ((double) hdp->samples_taken);
    int grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;

    int num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    double* distr;
    for (int id = 0; id < num_dps; id++){
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }

        distr = dp->posterior_predictive;

        for (int i = 0; i < grid_length; i++) {
            distr[i] = distr[i] * inv_sample_size;
        }

        dp->spline_slopes = spline_knot_slopes(grid, distr, grid_length);
    }
    
    hdp->splines_finalized = true;
}

double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int dp_id) {
    if (!hdp->splines_finalized) {
        fprintf(stderr, "Must finalize distributions before querying densities.\n");
        exit(EXIT_FAILURE);
    }

    if (dp_id < 0 || dp_id >= hdp->num_dps) {
        fprintf(stderr, "Hierarchical Dirichlet process has no Dirichlet process with this ID.\n");
        exit(EXIT_FAILURE);
    }

    DirichletProcess* dp = hdp->dps[dp_id];
    while (!dp->observed) {
        dp = dp->parent;
    }

    double interp =  grid_spline_interp(x, hdp->sampling_grid, dp->posterior_predictive,
                                        dp->spline_slopes, hdp->grid_length);
    if (interp > 0.0) {
        return interp;
    }
    else {
        return 0.0;
    }
}



//struct Factor {
//    FactorType factor_type;
//    struct Factor* parent;
//    stSet* children;
//    double* factor_data;
//    struct DirichletProcess* dp;
//};

void serialize_factor_tree_internal(FILE* out, Factor* fctr, int* parent_id_ptr, int* next_fctr_id, uintptr_t data_start) {
    int id = *next_fctr_id;
    (*next_fctr_id)++;
    if (parent_id_ptr == NULL) {
        fprintf(out, "-\t[");
    }
    else {
        fprintf(out, "%d\t[", *parent_id_ptr);
    }
    if (fctr->factor_type == DATA_PT) {
        uintptr_t data_pos = (uintptr_t) fctr->factor_data;
        fprintf(out, "%d", (data_pos - data_start) / sizeof(int));
    }
    else if (fctr->factor_type == BASE) {
        double* param_array = fctr->factor_data;
        for (int i = 0; i < N_IG_NUM_PARAMS; i++) {
            fprintf(out, "%15lf;", param_array[i]);
        }
        fprintf(out, "%15lf", param_array[N_IG_NUM_PARAMS]);
    }
    fprintf(out, "]\n");
    
    stSetIterator* iter = stSet_getIterator(fctr->children);
    Factor* child_fctr = (Factor*) stSet_getNext(iter);
    while (child_iter != NULL) {
        serialize_factor_tree_internal(out, child_fctr, &id, next_fctr_id, data_start);
        child_fctr = (Factor*) stSet_getNext(iter);
    }
    stSet_destructIterator(iter);
}

//struct HierarchicalDirichletProcess {
//    bool finalized;
//    double* data;
//    int* data_pt_dp_id;
//    int data_length;
//
//    struct DirichletProcess* base_dp;
//    struct DirichletProcess** dps;
//    int num_dps;
//
//    // normal-inverse gamma parameters
//    double mu;
//    double nu;
//    double two_alpha;
//    double beta;
//
//    double* sampling_grid;
//    int grid_length;
//
//    int samples_taken;
//    bool splines_finalized;
//
//    // TODO: replace this with my new offset log gamma memo
//    struct SumOfLogsMemo* log_sum_memo;
//
//    int depth;
//    bool sample_gamma;
//    double* gamma;
//
//    double* gamma_alpha;
//    double* gamma_beta;
//    double* w_aux_vector;
//    bool* s_aux_vector;
//};


void serialize_hdp(Hierarchical* hdp, const char* out_filepath) {
    FILE* out = fopen(out_filepath, "w");
    
    int num_data = hdp->num_data;
    double* data = hdp->data;
    int* dp_ids = hdp->data_pt_dp_id;
    int grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;
    int depth = hdp->depth;
    int gamma_params = hdp->gamma;
    int gamma_alpha = hdp->gamma_alpha;
    int gamma_beta = hdp->gamma_beta;
    double* w_aux_vector = hdp->w_aux_vector;
    bool* s_aux_vector = hdp->s_aux_vector;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* base_dp = hdp->base_dp;
    
    // finalized
    fprintf(out, "%d\n", (int) hdp->finalized);
    fprintf(out, "%d\n", (int) hdp->splines_finalized);
    // data
    for (int i = 0; i < num_data - 1; i++) {
        fprintf(out, "%15lf\t", data[i]);
    }
    fprintf(out, "%15lf\n", data[num_data - 1]);
    // dp ids
    for (int i = 0; i < num_data - 1; i++) {
        fprintf(out, "%d\t", dp_ids[i]);
    }
    fprintf(out, "%d\n", dp_ids[num_data - 1]);
    // base params
    fprintf(out, "%15lf\t%15lf\t%15lf\t%15lf\n", hdp->mu, hdp->nu, hdp->two_alpha, hdp->beta);
    // sampling grid
    fprintf(out, "%15lf\t%15l\t%d", grid[0], grid[grid_length - 1], grid_length);
    // gamma
    for (int i = 0; i < depth - 1; i++) {
        fprintf(out, "%15lf\t", gamma_params[i]);
    }
    fprintf(out, "%15lf\n", gamma_params[depth - 1]);
    // dirichlet processes
    DirichletProcess* dp;
    double* post_pred;
    double* slopes;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];
        //fprintf("%d\t", dp->num_factor_children);
        //fprintf("%d\t", (int) dp->observed);
        // children: no
        // gamma: no
        // depth: no
        // factors: no
        // num_children: no
        // wt: no
        // cached vals: no
        
        // parent
        if (dp == base_dp) {
            fprintf("-\t");
        }
        else {
            fprintf(out, "%d\t", dp->parent->id);
        }        post_pred = dp->posterior_predictive;
        slopes = dp->spline_slopes;
        // post pred
        fprintf(out, "[");
        for (int i = 0; i < grid_length - 1; i++) {
            fprintf(out, "%15lf;", post_pred[i]);
        }
        fprintf(out, "%15lf\n", post_pred[grid_length - 1]);
        // spline slopes
        fprintf(out, "]\t[");
        if (hdp->splines_finalized) {
            for (int i = 0; i < grid_length - 1; i++) {
                fprintf(out, "%15lf;", slopes[i]);
            }
            fprintf(out, "%15lf\n", slopes[grid_length - 1]);
        }
        fprintf(out, "]\n");
    }
    fprintf(out, "*\n");
    // gamma distr params
    if (hdp->sample_gamma) {
        // alpha
        for (int i = 0; i < depth - 1; i++) {
            fprintf(out, "%15lf\t", gamma_alpha[i]);
        }
        fprintf(out, "%15lf\n", gamma_alpha[depth - 1]);
        // beta
        for (int i = 0; i < depth - 1; i++) {
            fprintf(out, "%15lf\t", gamma_beta[i]);
        }
        fprintf(out, "%15lf\n", gamma_beta[depth - 1]);
        // w
        for (int i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%15lf\t", w_aux_vector[i]);
        }
        fprintf(out, "%15lf\n", w_aux_vector[num_dps - 1]);
        // s
        for (int i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%d\t", s_aux_vector[i]);
        }
        fprintf(out, "%d\n", s_aux_vector[num_dps - 1]);
    }
    else {
        fprintf(out, "\n");
    }
    // factors
    fprintf("*\n");
    int next_fctr_id = 0;
    uintptr_t data_start = (uintptr_t) hdp->data;
    
    stSetIterator* iter = stSet_getIterator(base_dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(iter);
    while (fctr != NULL) {
        serialize_factor_tree_internal(out, fctr, NULL, &next_fctr_id, data_start);
        fctr = (Factor*) stSet_getNext(iter);
    }
    stSet_destructIterator(iter);
    
    fclose(out);
}

HierarchicalDirichletProcess* deserialize_hdp(const char* filepath) {
    FILE* in = fopen(filepath, "r");
    
    // finalized
    char* end;
    char* line = stFile_getLineFromFile(in);
    bool finalized = (bool) strtol(line, &end, 10);
    free(line);
    // splines finalized
    line = stFile_getLineFromFile(in);
    bool splines_finalized = (bool) strtol(line, &end, 10);
    free(line);
    // data
    line = stFile_getLineFromFile(in);
    stList* tokens = stString_split(line);
    int data_length = stList_length(tokens);
    double* data = (double*) malloc(sizeof(double) * data_length);
    for (int i = 0; i < data_length; i++) {
        sscanf(stList_get(tokens, i), "%lf", &(data[i]));
    }
    free(line);
    stList_destruct(tokens);
    // dp ids
    line = stFile_getLineFromFile(in);
    tokens = stString_split(line);
    int* dp_ids = (int*) malloc(sizeof(int) * data_length);
    for (int i = 0; i < data_length; i++) {
        sscanf(stList_get(tokens, i), "%d", &(dp_ids[i]));
    }
    free(line);
    stList_destruct(tokens);
    // base params
    line = stFile_getLineFromFile(in);
    double mu, nu, two_alpha, beta;
    sscanf(line, "%15lf\t%15lf\t%15lf\t%15lf", &mu, &nu, &two_alpha, &beta);
    free(line);
    // sampling grid
    line = stFile_getLineFromFile(in);
    double grid_start, grid_stop;
    int grid_length;
    sscanf(line, "%15lf\t%15l\t%d", &grid_start, &grid_stop, &grid_length);
    free(line);
    // gamma
    line = stFile_getLineFromFile(in);
    tokens = stString_split(line);
    int depth = stList_length(tokens);
    double* gamma_params = (double*) malloc(sizeof(double) * depth);
    for (int i = 0; i < depth; i++) {
        sscanf(stList_get(tokens, i), "%lf", &(gamma_params[i]));
    }
    free(line);
    stList_destruct(tokens);
    // gamma distr params
    line = stFile_getLineFromFile(in);
    tokens = stString_split(line);
    bool sample_gamma = false;
    double* gamma_alpha;
    double* gamma_beta;
    double* w;
    bool* s;
    int s_int;
    int num_dps; // TODO: careful with this, don't like it
    if (stList_length(tokens) > 0) {
        // gamma alpha
        sample_gamma = true;
        gamma_alpha = (double*) malloc(sizeof(double) * depth);
        for (int i = 0; i < depth; i++) {
            sscanf(stList_get(tokens, i), "%lf", &(gamma_alpha[i]));
        }
        free(line);
        stList_destruct(tokens);
        // gamma beta
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        gamma_beta = (double*) malloc(sizeof(double) * depth);
        for (int i = 0; i < depth; i++) {
            sscanf(stList_get(tokens, i), "%lf", &(gamma_beta[i]));
        }
        free(line);
        stList_destruct(tokens);
        // w
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        num_dps = stList_length(tokens);
        w = (int*) malloc(sizeof(int) * num_dps);
        for (int i = 0; i < num_dps; i++) {
            sscanf(stList_get(tokens, i), "%lf", &(w[i]));
        }
        free(line);
        stList_destruct(tokens);
        // s
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        s = (bool*) malloc(sizeof(bool) * num_dps);
        for (int i = 0; i < num_dps; i++) {
            sscanf(stList_get(tokens, i), "%d", &s_int);
            s[i] = (bool) s_int;
        }
    }
    free(line);
    stList_destruct(tokens);
    
    
    
    //TODO: finish here
    
    
    
    
    
    
    // dirichlet processes
    DirichletProcess* dp;
    double* post_pred;
    double* slopes;
    for (int i = 0; i < num_dps; i++) {
        dp = dps[i];
        //fprintf("%d\t", dp->num_factor_children);
        //fprintf("%d\t", (int) dp->observed);
        // children: no
        // gamma: no
        // depth: no
        // factors: no
        // num_children: no
        // wt: no
        // cached vals: no
        
        // parent
        if (dp == base_dp) {
            fprintf("-\t");
        }
        else {
            fprintf(out, "%d\t", dp->parent->id);
        }        post_pred = dp->posterior_predictive;
        slopes = dp->spline_slopes;
        // post pred
        fprintf(out, "[");
        for (int i = 0; i < grid_length - 1; i++) {
            fprintf(out, "%15lf;", post_pred[i]);
        }
        fprintf(out, "%15lf\n", post_pred[grid_length - 1]);
        // spline slopes
        fprintf(out, "]\t[");
        if (hdp->splines_finalized) {
            for (int i = 0; i < grid_length - 1; i++) {
                fprintf(out, "%15lf;", slopes[i]);
            }
            fprintf(out, "%15lf\n", slopes[grid_length - 1]);
        }
        fprintf(out, "]\n");
    }
    // gamma distr params
    if (hdp->sample_gamma) {
        // alpha
        for (int i = 0; i < depth - 1; i++) {
            fprintf(out, "%15lf\t", gamma_alpha[i]);
        }
        fprintf(out, "%15lf\n", gamma_alpha[depth - 1]);
        // beta
        for (int i = 0; i < depth - 1; i++) {
            fprintf(out, "%15lf\t", gamma_beta[i]);
        }
        fprintf(out, "%15lf\n", gamma_beta[depth - 1]);
        // w
        for (int i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%15lf\t", w_aux_vector[i]);
        }
        fprintf(out, "%15lf\n", w_aux_vector[num_dps - 1]);
        // s
        for (int i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%d\t", s_aux_vector[i]);
        }
        fprintf(out, "%d\n", s_aux_vector[num_dps - 1]);
    }
    else {
        fprintf(out, "\n");
    }
    // factors
    fprintf("*\n");
    int next_fctr_id = 0;
    uintptr_t data_start = (uintptr_t) hdp->data;
    
    stSetIterator* iter = stSet_getIterator(base_dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(iter);
    while (fctr != NULL) {
        serialize_factor_tree_internal(out, fctr, NULL, &next_fctr_id, data_start);
        fctr = (Factor*) stSet_getNext(iter);
    }
    stSet_destructIterator(iter);
    
    fclose(out);

}

//struct DirichletProcess {
//    int id;
//    struct HierarchicalDirichletProcess* hdp;
//    double* gamma;
//    int depth;
//
//    struct DirichletProcess* parent;
//    stList* children;
//    stSet* factors;
//    int num_factor_children;
//
//    double base_factor_wt;
//    double* posterior_predictive;
//    double* spline_slopes;
//
//    double cached_factor_mean;
//    double cached_factor_sum_sq_dev;
//    int cached_factor_size;
//
//    bool observed;
//};












