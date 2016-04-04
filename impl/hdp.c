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

#ifndef MINUS_INF
#define MINUS_INF -0.5 * DBL_MAX
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338
#endif

typedef struct DirichletProcess {
    int64_t id;
    struct HierarchicalDirichletProcess* hdp;
    double* gamma;
    int64_t depth;
    
    struct DirichletProcess* parent;
    stList* children;
    stSet* factors;
    int64_t num_factor_children;
    
    double base_factor_wt;
    double* posterior_predictive;
    double* spline_slopes;
    
    // TODO: do something about these cached variables. maybe one store for the entire hdp?
    double cached_factor_mean;
    double cached_factor_sum_sq_dev;
    int64_t cached_factor_size;
    
    bool observed;
} DirichletProcess;

typedef enum FactorType {
    BASE,
    MIDDLE,
    DATA_PT,
    FUZZY_DATA_PT
} FactorType;

typedef struct Factor {
    FactorType factor_type;
    void* factor_data;
} Factor;

typedef struct BaseFactorData {
    stSet* children;
    struct DirichletProcess* dp;
    double mu;
    double nu;
    double alpha;
    double beta;
    double log_posterior_term;
} BaseFactorData;

typedef struct MiddleFactorData {
    struct Factor* parent;
    stSet* children;
    struct DirichletProcess* dp;
} MiddleFactorData;

typedef struct DataPtFactorData {
    struct Factor* parent;
    double data_pt;
} DataPtFactorData;

typedef struct FuzzyDataPtFactorData {
    struct Factor* parent;
    double data_pt;
    DirichletProcess** fuzzy_dps;
    double* fuzzy_dp_cdf;
    int64_t num_fuzzy_dps;
} FuzzyDataPtFactorData;


// TEMPLATE:
//switch (fctr->factor_type) {
//    case BASE:
//    {
//        BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
//        break;
//    }
//    case MIDDLE:
//    {
//        MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
//        break;
//    }
//    case DATA_PT:
//    {
//        DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
//        break;
//    }
//    case FUZZY_DATA_PT:
//    {
//        FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
//        break;
//    }
//    default:
//    {
//        fprintf(stderr, "Unsupported factor type.\n");
//        exit(EXIT_FAILURE);
//    }
//}


struct HierarchicalDirichletProcess {
    bool finalized;
    bool fuzzy_assignments;
    bool has_data;
    
    double* data;
    int64_t data_length;
    
    int64_t* data_pt_dp_id;
    
    int64_t** data_pt_fuzzy_dp_ids;
    double** data_pt_fuzzy_dp_probs;
    int64_t* data_pt_num_fuzzy_dps;

    struct DirichletProcess* base_dp;
    struct DirichletProcess** dps;
    int64_t num_dps;

    // normal-inverse gamma parameters
    double mu;
    double nu;
    double alpha;
    double beta;

    double* sampling_grid;
    int64_t grid_length;
    
    int64_t samples_taken;
    bool splines_finalized;
    
    int64_t depth;
    bool sample_gamma;
    double* gamma;
    
    double* gamma_alpha;
    double* gamma_beta;
    double* w_aux_vector;
    bool* s_aux_vector;
    
    stSet* distr_metric_memos;
};

struct DistributionMetricMemo {
    int64_t num_distrs;
    double* memo_matrix;
    
    HierarchicalDirichletProcess* hdp;
    double (*metric_func) (HierarchicalDirichletProcess*, int64_t, int64_t);
};

bool is_structure_finalized(HierarchicalDirichletProcess* hdp) {
    return hdp->finalized;
}

bool is_gamma_random(HierarchicalDirichletProcess* hdp) {
    return hdp->sample_gamma;
}

bool is_sampling_finalized(HierarchicalDirichletProcess* hdp) {
    return hdp->splines_finalized;
}

bool are_data_pt_assignments_fuzzy(HierarchicalDirichletProcess* hdp) {
    return hdp->fuzzy_assignments;
}

int64_t get_num_dir_proc(HierarchicalDirichletProcess* hdp) {
    return hdp->num_dps;
}

int64_t get_depth(HierarchicalDirichletProcess* hdp) {
    return hdp->depth;
}

int64_t get_num_data(HierarchicalDirichletProcess* hdp) {
    return hdp->data_length;
}

//double* get_data_copy(HierarchicalDirichletProcess* hdp) {
//    int64_t data_length = hdp->data_length;
//    double* data = (double*) malloc(sizeof(double) * data_length);
//    for (int64_t i = 0; i < data_length; i++) {
//        data[i] = hdp->data[i];
//    }
//    return data;
//}
//
//int64_t* get_data_pt_dp_ids_copy(HierarchicalDirichletProcess* hdp) {
//    int64_t data_length = hdp->data_length;
//    int64_t* dp_ids = (int64_t*) malloc(sizeof(int64_t) * data_length);
//    for (int64_t i = 0; i < data_length; i++) {
//        dp_ids[i] = hdp->data_pt_dp_id[i];
//    }
//    return dp_ids;
//}

double* get_gamma_params_copy(HierarchicalDirichletProcess* hdp) {
    int64_t depth = hdp->depth;
    double* gamma_params = (double*) malloc(sizeof(double) * depth);
    for (int64_t i = 0; i < depth; i++) {
        gamma_params[i] = hdp->gamma[i];
    }
    return gamma_params;
}

double get_mu(HierarchicalDirichletProcess* hdp) {
    return hdp->mu;
}

double get_nu(HierarchicalDirichletProcess* hdp) {
    return hdp->nu;
}

double get_alpha(HierarchicalDirichletProcess* hdp) {
    return hdp->alpha;
}

double get_beta(HierarchicalDirichletProcess* hdp) {
    return hdp->beta;
}

int64_t get_grid_length(HierarchicalDirichletProcess* hdp) {
    return hdp->grid_length;
}

double* get_sampling_grid_copy(HierarchicalDirichletProcess* hdp) {
    int64_t grid_length = hdp->grid_length;
    double* sampling_grid = (double*) malloc(sizeof(double) * grid_length);
    for (int64_t i = 0; i < grid_length; i++) {
        sampling_grid[i] = hdp->sampling_grid[i];
    }
    return sampling_grid;
}

double* get_gamma_alpha_params_copy(HierarchicalDirichletProcess* hdp) {
    if (!hdp->sample_gamma) {
        fprintf(stderr, "Hierarchical Dirichlet process is not sampling gamma parameters.");
        exit(EXIT_FAILURE);
    }
    int64_t depth = hdp->depth;
    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
    for (int64_t i = 0; i < depth; i++) {
        gamma_alpha[i] = hdp->gamma_alpha[i];
    }
    return gamma_alpha;
}

double* get_gamma_beta_params_copy(HierarchicalDirichletProcess* hdp) {
    if (!hdp->sample_gamma) {
        fprintf(stderr, "Hierarchical Dirichlet process is not sampling gamma parameters.");
        exit(EXIT_FAILURE);
    }
    int64_t depth = hdp->depth;
    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
    for (int64_t i = 0; i < depth; i++) {
        gamma_beta[i] = hdp->gamma_beta[i];
    }
    return gamma_beta;
}

int64_t get_dir_proc_num_factors(HierarchicalDirichletProcess* hdp, int64_t dp_id) {
    if (dp_id < 0 || dp_id >= hdp->num_dps) {
        fprintf(stderr, "Hierarchical Dirichlet process has no Dirichlet process with this ID.\n");
        exit(EXIT_FAILURE);
    }
    
    DirichletProcess* dp = hdp->dps[dp_id];
    return stSet_size(dp->factors);
}

int64_t get_dir_proc_parent_id(HierarchicalDirichletProcess* hdp, int64_t dp_id) {
    if (dp_id < 0 || dp_id >= hdp->num_dps) {
        fprintf(stderr, "Hierarchical Dirichlet process has no Dirichlet process with this ID.\n");
        exit(EXIT_FAILURE);
    }
    
    DirichletProcess* dp = hdp->dps[dp_id];
    if (dp->parent == NULL) {
        return -1;
    }
    else {
        return dp->parent->id;
    }
}

void get_dp_depths_internal(int64_t* dp_depths, DirichletProcess* dp, int64_t depth) {
    dp_depths[dp->id] = depth;
    
    stListIterator* dp_child_iter = stList_getIterator(dp->children);
    DirichletProcess* dp_child = stList_getNext(dp_child_iter);
    while (dp_child != NULL) {
        get_dp_depths_internal(dp_depths, dp_child, depth + 1);
        dp_child = stList_getNext(dp_child_iter);
    }
    stList_destructIterator(dp_child_iter);
}

int64_t* get_dp_depths(HierarchicalDirichletProcess* hdp) {
    int64_t* dp_depths = (int64_t*) malloc(sizeof(int64_t) * hdp->num_dps);
    get_dp_depths_internal(dp_depths, hdp->base_dp, 0);
    return dp_depths;
}

DistributionMetricMemo* new_distr_metric_memo(HierarchicalDirichletProcess* hdp,
                                              double (*metric_func) (HierarchicalDirichletProcess*, int64_t, int64_t)) {
    DistributionMetricMemo* memo = (DistributionMetricMemo*) malloc(sizeof(DistributionMetricMemo));
    int64_t num_dps = hdp->num_dps;
    memo->num_distrs = num_dps;
    
    int64_t num_entries = ((num_dps - 1) * num_dps) / 2;
    double* memo_matrix = (double*) malloc(sizeof(double) * num_entries);
    memo->memo_matrix = memo_matrix;
    for (int64_t i = 0; i < num_entries; i++) {
        memo_matrix[i] = -1.0;
    }
    
    memo->hdp = hdp;
    memo->metric_func = metric_func;
    
    stSet_insert(hdp->distr_metric_memos, memo);
    
    return memo;
}

void destroy_distr_metric_memo(void* memo) {
    DistributionMetricMemo* metric_memo = (DistributionMetricMemo*) memo;
    free(metric_memo->memo_matrix);
    free(metric_memo);
}

void cache_base_factor_params(Factor* fctr, double mu, double nu, double alpha, double beta, double log_post_term) {
    if (fctr->factor_type != BASE) {
        fprintf(stderr, "Can only cache parameters for base factors.\n");
        exit(EXIT_FAILURE);
    }
    
    BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;

    fctr_data->mu = mu;
    fctr_data->nu = nu;
    fctr_data->alpha = alpha;
    fctr_data->beta = beta;
    fctr_data->log_posterior_term = log_post_term;
}

Factor* new_base_factor(HierarchicalDirichletProcess* hdp) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = BASE;
    
    BaseFactorData* fctr_data = (BaseFactorData*) malloc(sizeof(BaseFactorData));
    fctr_data->children = stSet_construct();
    fctr_data->dp = hdp->base_dp;
    fctr->factor_data = (void*) fctr_data;
    cache_base_factor_params(fctr, hdp->mu, hdp->nu, hdp->alpha, hdp->beta, 1.0);

    stSet_insert(hdp->base_dp->factors, (void*) fctr);

    return fctr;
}

Factor* new_middle_factor(DirichletProcess* dp) {
    if (dp->parent == NULL) {
        fprintf(stderr, "Attempted to create middle factor in root Dirichlet process.\n");
        exit(EXIT_FAILURE);
    }

    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = MIDDLE;
    
    MiddleFactorData* fctr_data = (MiddleFactorData*) malloc(sizeof(MiddleFactorData));
    fctr_data->dp = dp;
    fctr_data->children =stSet_construct();
    fctr_data->parent = NULL; // note: assigning to parent handled externally
    
    fctr->factor_data = (void*) fctr_data;

    stSet_insert(dp->factors, (void*) fctr);
    
    return fctr;
}

Factor* new_data_pt_factor(double data_pt) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = DATA_PT;
    
    DataPtFactorData* fctr_data = (DataPtFactorData*) malloc(sizeof(DataPtFactorData));
    fctr_data->parent = NULL; // note: assigning to parent handled externally
    fctr_data->data_pt = data_pt;
    
    fctr->factor_data = (void*) fctr_data;

    return fctr;
}

Factor* new_fuzzy_data_pt_factor(double data_pt, DirichletProcess** fuzzy_dps, double* fuzzy_dp_probs,
                                 int64_t num_fuzzy_dps) {
    Factor* fctr = (Factor*) malloc(sizeof(Factor));
    fctr->factor_type = FUZZY_DATA_PT;
    
    FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) malloc(sizeof(FuzzyDataPtFactorData));
    fctr_data->parent = NULL; // note: assigning to parent handled externally
    fctr_data->data_pt = data_pt;
    fctr_data->fuzzy_dps = fuzzy_dps;
    fctr_data->fuzzy_dp_cdf = fuzzy_dp_probs;
    fctr_data->num_fuzzy_dps = num_fuzzy_dps;
    
    for (int64_t i = 1; i < num_fuzzy_dps; i++) {
        fuzzy_dp_probs[i] = fuzzy_dp_probs[i] + fuzzy_dp_probs[i - 1];
    }
    
    fctr->factor_data = (void*) fctr_data;
    
    return fctr;
}

DirichletProcess* get_factor_dir_proc(Factor* fctr) {
    DirichletProcess* dp;
    switch (fctr->factor_type) {
        case BASE:
        {
            BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
            dp = fctr_data->dp;
            break;
        }
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            dp = fctr_data->dp;
            break;
        }
        default:
        {
            fprintf(stderr, "Attempted to access Dirichlet process of leaf factor.\n");
            exit(EXIT_FAILURE);
        }
    }
    return dp;
}

stSet* get_factor_children(Factor* fctr) {
    stSet* children;
    switch (fctr->factor_type) {
        case BASE:
        {
            BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
            children = fctr_data->children;
            break;
        }
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            children = fctr_data->children;
            break;
        }
        default:
        {
            fprintf(stderr, "Attempted to access children of leaf factor.\n");
            exit(EXIT_FAILURE);
        }
    }
    return children;
}

void destroy_factor(Factor* fctr) {
    
    Factor* parent = NULL;
    stSet* children = NULL;
    DirichletProcess* dp = NULL;
    //printf("finding factor type:\n");
    switch (fctr->factor_type) {
        case BASE:
        {
            //printf("base\n");
            BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
            dp = fctr_data->dp;
            children = fctr_data->children;
            break;
        }
            
        case MIDDLE:
        {
            //printf("middle\n");
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            dp = fctr_data->dp;
            parent = fctr_data->parent;
            children = fctr_data->children;
            break;
        }
            
        case DATA_PT:
        {
            //printf("data pt\n");
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            break;
        }
        
        case FUZZY_DATA_PT:
        {
            //printf("fuzzy\n");
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            free(fctr_data->fuzzy_dps);
            free(fctr_data->fuzzy_dp_cdf);
            break;
        }
            
        default:
        {
            fprintf(stderr, "Attempted to destroy unsupported factor type.\n");
            exit(EXIT_FAILURE);
            break;
        }
    }
    
    //printf("free factor data\n");
    free(fctr->factor_data);
    
    if (children != NULL) {
        //printf("destroy child set\n");
        if (stSet_size(children) > 0) {
            fprintf(stderr, "Attempted to destroy factor that still has children.\n");
            exit(EXIT_FAILURE);
        }
        stSet_destruct(children);
    }

    if (parent != NULL) {
        //printf("updating parent\n");
        stSet* parents_children = get_factor_children(parent);
        //printf("removing from parent set\n");
        stSet_remove(parents_children, (void*) fctr);
        DirichletProcess* parent_dp = get_factor_dir_proc(parent);
        //printf("decrementing\n");
        (parent_dp  ->num_factor_children)--;
        if (stSet_size(parents_children) == 0) {
            destroy_factor(parent);
        }
    }

    if (dp != NULL) {
        //printf("removing from dp\n");
        stSet_remove(dp->factors, (void*) fctr);
    }

    free(fctr);
}

DirichletProcess* sample_fuzzy_data_pt_dp_assignment(Factor* fctr) {
    if (fctr->factor_type != FUZZY_DATA_PT) {
        fprintf(stderr, "Can only sample Dirichlet process assignment of fuzzy data point factors.\n");
        exit(EXIT_FAILURE);
    }
    
    FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
    int64_t num_fuzzy_dps = fctr_data->num_fuzzy_dps;
    DirichletProcess** fuzzy_dps = fctr_data->fuzzy_dps;
    double* cdf = fctr_data->fuzzy_dp_cdf;
    double total_prob = cdf[num_fuzzy_dps - 1];
    
    double draw = rand_uniform(total_prob);
    DirichletProcess* dp_choice = fuzzy_dps[bisect_left(draw, cdf, num_fuzzy_dps)];
    
    return dp_choice;
}

Factor* get_base_factor(Factor* fctr) {
    while (true) {
        switch (fctr->factor_type) {
            case BASE:
            {
                return fctr;
            }
            case MIDDLE:
            {
                MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
                fctr = fctr_data->parent;
                break;
            }
            case DATA_PT:
            {
                DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
                fctr = fctr_data->parent;
                break;
            }
            case FUZZY_DATA_PT:
            {
                FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
                fctr = fctr_data->parent;
                break;
            }
            default:
            {
                fprintf(stderr, "Unsupported factor type.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    return NULL;
}

double get_factor_data_pt(Factor* fctr) {
    double data_pt;
    switch (fctr->factor_type) {
        case DATA_PT:
        {
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            data_pt = fctr_data->data_pt;
            break;
        }
        case FUZZY_DATA_PT:
        {
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            data_pt = fctr_data->data_pt;
            break;
        }
        default:
        {
            fprintf(stderr, "Attempted to access data point from non-leaf factor.\n");
            exit(EXIT_FAILURE);
        }
    }

    return data_pt;
}

Factor* get_factor_parent(Factor* fctr) {
    Factor* parent;
    switch (fctr->factor_type) {
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            break;
        }
        case DATA_PT:
        {
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            break;
        }
        case FUZZY_DATA_PT:
        {
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            break;
        }
        default:
        {
            fprintf(stderr, "Attempted to access factor parent of base factor.\n");
            exit(EXIT_FAILURE);
        }
    }
    return parent;
}

Factor** get_factor_parent_ptr(Factor* fctr) {
    Factor** parent_ptr;
    switch (fctr->factor_type) {
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            parent_ptr = &(fctr_data->parent);
            break;
        }
        case DATA_PT:
        {
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            parent_ptr = &(fctr_data->parent);
            break;
        }
        case FUZZY_DATA_PT:
        {
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            parent_ptr = &(fctr_data->parent);
            break;
        }
        default:
        {
            fprintf(stderr, "Attempted to access factor parent of base factor.\n");
            exit(EXIT_FAILURE);
        }
    }
    return parent_ptr;
}

void get_factor_sum_internal(Factor* fctr, double* sum, int64_t* num_data) {
    if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
        *sum += get_factor_data_pt(fctr);
        // TODO: there should be a way to use the parent's counters instead of recounting the data pts
        (*num_data)++;
    }
    else {
        stSetIterator* child_iter = stSet_getIterator(get_factor_children(fctr));
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            get_factor_sum_internal(child_fctr, sum, num_data);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        stSet_destructIterator(child_iter);
    }
}

void get_factor_ssd_internal(Factor* fctr, double center, double* sum_sq_dev) {
    if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
        double dev = get_factor_data_pt(fctr) - center;
        *sum_sq_dev += dev * dev;
    }
    else {
        stSetIterator* child_iter = stSet_getIterator(get_factor_children(fctr));
        Factor* child_fctr = (Factor*) stSet_getNext(child_iter);
        while (child_fctr != NULL) {
            get_factor_ssd_internal(child_fctr, center, sum_sq_dev);
            child_fctr = (Factor*) stSet_getNext(child_iter);
        }
        stSet_destructIterator(child_iter);
    }
}

void get_factor_stats(Factor* fctr, double* mean_out, double* sum_sq_dev_out, int64_t* num_data_out) {
    *mean_out = 0.0;
    *sum_sq_dev_out = 0.0;
    *num_data_out = 0;
    get_factor_sum_internal(fctr, mean_out, num_data_out);
    *mean_out /= (double) *num_data_out;
    get_factor_ssd_internal(fctr, *mean_out, sum_sq_dev_out);
}

void add_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    if (fctr->factor_type != BASE) {
        fprintf(stderr, "Attempted to update posterior parameters for non-base factor.\n");
        exit(EXIT_FAILURE);
    }
    
    BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;

    double mu_prev = fctr_data->mu;
    double nu_prev = fctr_data->nu;
    double alpha_prev = fctr_data->alpha;
    double beta_prev = fctr_data->beta;

    double nu_post = nu_prev + num_data;
    double mu_post = (mu_prev * nu_prev + mean * num_data) / nu_post;
    double alpha_post = alpha_prev + 0.5 * num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_post = beta_prev + .5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_post, alpha_post, beta_post);

    cache_base_factor_params(fctr, mu_post, nu_post, alpha_post, beta_post, log_post_term);
}

void remove_update_base_factor_params(Factor* fctr, double mean, double sum_sq_devs, double num_data) {
    if (fctr->factor_type != BASE) {
        fprintf(stderr, "Attempted to update posterior parameters for non-base factor.\n");
        exit(EXIT_FAILURE);
    }
    BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;

    double mu_post = fctr_data->mu;
    double nu_post = fctr_data->nu;
    double alpha_post = fctr_data->alpha;
    double beta_post = fctr_data->beta;

    double nu_prev = nu_post - num_data;
    double mu_prev = (mu_post * nu_post - mean * num_data) / nu_prev;
    double alpha_prev = alpha_post - 0.5 * num_data;

    double mean_dev = mean - mu_prev;
    double sq_mean_dev = nu_prev * num_data * mean_dev * mean_dev / nu_post;

    double beta_prev = beta_post - 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_post_term = log_posterior_conditional_term(nu_prev, alpha_prev, beta_prev);

    cache_base_factor_params(fctr, mu_prev, nu_prev, alpha_prev, beta_prev, log_post_term);
}

double factor_parent_joint_log_likelihood(Factor* fctr, Factor* parent) {
    Factor* base_fctr = get_base_factor(parent);
    DirichletProcess* dp = get_factor_dir_proc(fctr);

    double num_reassign = (double) dp->cached_factor_size;
    double mean_reassign = dp->cached_factor_mean;
    double sum_sq_devs = dp->cached_factor_sum_sq_dev;

    BaseFactorData* base_fctr_data = (BaseFactorData*) base_fctr->factor_data;
    
    double mu_denom = base_fctr_data->mu;
    double nu_denom = base_fctr_data->nu;
    double alpha_denom = base_fctr_data->alpha;
    double beta_denom = base_fctr_data->beta;

    double nu_numer = nu_denom + num_reassign;
    double alpha_numer = alpha_denom + 0.5 * num_reassign;

    double mean_dev = mean_reassign - mu_denom;
    double sq_mean_dev = nu_denom * num_reassign * mean_dev * mean_dev / nu_numer;

    double beta_numer = beta_denom + 0.5 * (sum_sq_devs + sq_mean_dev);

    double log_denom = base_fctr_data->log_posterior_term;
    double log_numer = log_posterior_conditional_term(nu_numer, alpha_numer, beta_numer);

    return -0.5 * num_reassign * log(2.0 * M_PI) + log_numer - log_denom;
}

double data_pt_factor_parent_likelihood(Factor* data_pt_fctr, Factor* parent) {
    if (data_pt_fctr->factor_type != DATA_PT && data_pt_fctr->factor_type != FUZZY_DATA_PT) {
        fprintf(stderr, "Can only access data point likelihood for data point factors.\n");
        exit(EXIT_FAILURE);
    }
    
    double data_pt = get_factor_data_pt(data_pt_fctr);
    Factor* base_fctr = get_base_factor(parent);
    
    BaseFactorData* base_fctr_data = (BaseFactorData*) base_fctr->factor_data;
    
    double mu_denom = base_fctr_data->mu;
    double nu_denom = base_fctr_data->nu;
    double alpha_denom = base_fctr_data->alpha;
    double beta_denom = base_fctr_data->beta;

    double nu_numer = nu_denom + 1.0;

    double mean_dev = data_pt - mu_denom;
    double sq_mean_dev = nu_denom * mean_dev * mean_dev / nu_numer;

    double alpha_numer = alpha_denom + 0.5;
    double beta_numer = beta_denom + 0.5 * sq_mean_dev;
    
    double log_denom = base_fctr_data->log_posterior_term;
    double log_numer = log_posterior_conditional_term(nu_numer, alpha_numer, beta_numer);

    return (1.0 / sqrt(2.0 * M_PI)) * exp(log_numer - log_denom);
}

void evaluate_posterior_predictive(Factor* base_fctr, double* x, double* pdf_out, int64_t length) {
    if (base_fctr->factor_type != BASE) {
        fprintf(stderr, "Can only evaluate posterior predictive of base factors.\n");
        exit(EXIT_FAILURE);
    }
    
    BaseFactorData* base_fctr_data = (BaseFactorData*) base_fctr->factor_data;
    
    double mu_denom = base_fctr_data->mu;
    double nu_denom = base_fctr_data->nu;
    double alpha_denom = base_fctr_data->alpha;
    double beta_denom = base_fctr_data->beta;
    
    double log_denom = base_fctr_data->log_posterior_term;

    double nu_numer = nu_denom + 1.0;
    double alpha_numer = alpha_denom + 0.5;
    double nu_ratio = nu_denom / nu_numer;
    double pi_factor = 1.0 / sqrt(2.0 * M_PI);

    double mean_dev;
    double sq_mean_dev;
    double beta_numer;
    double log_numer;
    
    for (int64_t i = 0; i < length; i++) {
        mean_dev = x[i] - mu_denom;
        sq_mean_dev = nu_ratio * mean_dev * mean_dev;
        beta_numer = beta_denom + 0.5 * sq_mean_dev;

        log_numer = log_posterior_conditional_term(nu_numer, alpha_numer, beta_numer);

        pdf_out[i] = pi_factor * exp(log_numer - log_denom);
    }
}

void evaluate_prior_predictive(HierarchicalDirichletProcess* hdp,
                               double* x, double* pdf_out, int64_t length) {
    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    double alpha = hdp->alpha;
    double beta = hdp->beta;

    double nu_factor = nu / (2.0 * (nu + 1.0) * beta);
    double alpha_term = exp(lgamma(alpha + 0.5) - lgamma(alpha));
    double beta_term = sqrt(nu_factor / M_PI);
    double constant_term = alpha_term * beta_term;
    double alpha_power = -alpha - 0.5;
    
    for (int64_t i = 0; i < length; i++) {
        double dev = x[i] - mu;
        double var_term = pow(1.0 + nu_factor * dev * dev, alpha_power);

        pdf_out[i] = constant_term * var_term;
    }
}

double prior_likelihood(HierarchicalDirichletProcess* hdp, Factor* fctr) {
    if (fctr->factor_type != DATA_PT && fctr->factor_type != FUZZY_DATA_PT) {
        fprintf(stderr, "Cannot calculate point prior likelihood from non-data point factor.\n");
    }

    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    double alpha = hdp->alpha;
    double beta = hdp->beta;

    double data_pt = get_factor_data_pt(fctr);
    double dev = data_pt - mu;

    double alpha_term = exp(lgamma(alpha + 0.5) - lgamma(alpha));
    double nu_term = nu / (2.0 * (nu + 1.0) * beta);
    double beta_term = pow(1.0 + nu_term * dev * dev, -alpha - 0.5);
    
    return alpha_term * sqrt(nu_term / M_PI) * beta_term;
}

double prior_joint_log_likelihood(HierarchicalDirichletProcess* hdp, Factor* fctr) {
    if (fctr->factor_type != MIDDLE) {
        fprintf(stderr, "Cannot calculate joint prior likelihood from non-middle factor.\n");
    }
    
    //TODO: this could be made more efficient with some precomputed variables stashed in HDP
    double mu = hdp->mu;
    double nu = hdp->nu;
    double alpha = hdp->alpha;
    double beta = hdp->beta;
    
    DirichletProcess* dp = get_factor_dir_proc(fctr);
    double num_reassign = (double) dp->cached_factor_size;;
    double mean_reassign = dp->cached_factor_mean;
    double sum_sq_devs = dp->cached_factor_sum_sq_dev;
    
    double mean_dev = mean_reassign - mu;
    double sq_mean_dev = nu * num_reassign * mean_dev * mean_dev / (nu + num_reassign);
    
    double log_alpha_term = lgamma(alpha + .5 * num_reassign) - lgamma(alpha);
    double log_nu_term = 0.5 * (log(nu) - log(nu + num_reassign));
    double log_pi_term = 0.5 * num_reassign * log(2.0 * M_PI);
    double log_beta_term_1 = alpha * log(beta);
    double log_beta_term_2 = (alpha + 0.5 * num_reassign)
                              * log(beta + 0.5 * (sum_sq_devs + sq_mean_dev));
    return log_alpha_term + log_nu_term - log_pi_term + log_beta_term_1 - log_beta_term_2;
}

// TODO: figure out how to break into chunks and spin up threads to reduce the sum behind the iterator
double unobserved_factor_likelihood(Factor* fctr, DirichletProcess* dp) {
    DirichletProcess* parent_dp = dp->parent;
    if (parent_dp == NULL) {
        return prior_likelihood(dp->hdp, fctr);
    }
    else {
        double parent_gamma = *(parent_dp->gamma);
        double likelihood = 0.0;
        double next_height_unobs_likelihood;
        int64_t num_parent_fctrs = stSet_size(parent_dp->factors);
        Factor** parent_fctrs = (Factor**) malloc(sizeof(Factor*) * num_parent_fctrs);
        
        #pragma omp parallel shared(likelihood,next_height_unobs_likelihood,parent_dp,num_parent_fctrs,parent_fctrs)
        {
            #pragma omp single nowait
            next_height_unobs_likelihood = unobserved_factor_likelihood(fctr, parent_dp);
            
             #pragma omp single
            {
                stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
                for (int64_t i = 0; i < num_parent_fctrs; i++) {
                    parent_fctrs[i] = (Factor*) stSet_getNext(parent_fctr_iter);
                }
                stSet_destructIterator(parent_fctr_iter);
            }
            
            double local_likelihood = 0.0;
            Factor* parent_fctr;
            #pragma omp for nowait
            for (int64_t i = 0; i < num_parent_fctrs; i++) {
                parent_fctr = parent_fctrs[i];
                local_likelihood += stSet_size(get_factor_children(parent_fctr))
                                    * data_pt_factor_parent_likelihood(fctr, parent_fctr);
            }
            
            #pragma omp critical
            likelihood += local_likelihood;
        }
        free(parent_fctrs);
        
        likelihood += parent_gamma * next_height_unobs_likelihood;
        
        likelihood /= (parent_gamma + (double) parent_dp->num_factor_children);
        
        return likelihood;
    }
}

//double unobserved_factor_likelihood(Factor* fctr, DirichletProcess* dp) {
//    DirichletProcess* parent_dp = dp->parent;
//    if (parent_dp == NULL) {
//        return prior_likelihood(dp->hdp, fctr);
//    }
//    else {
//        double parent_gamma = *(parent_dp->gamma);
//        double likelihood = 0.0;
//        
//        stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
//        
//        Factor* parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
//        double fctr_size;
//        while (parent_fctr != NULL) {
//            fctr_size = (double) stSet_size(parent_fctr->children);
//            likelihood += fctr_size * data_pt_factor_parent_likelihood(fctr, parent_fctr);
//            parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
//        }
//        stSet_destructIterator(parent_fctr_iter);
//        
//        likelihood += parent_gamma * unobserved_factor_likelihood(fctr, parent_dp);
//        
//        likelihood /= (parent_gamma + (double) parent_dp->num_factor_children);
//        
//        return likelihood;
//    }
//}


double unobserved_factor_joint_log_likelihood(Factor* fctr, DirichletProcess* dp) {
    DirichletProcess* parent_dp = dp->parent;
    if (parent_dp == NULL) {
        return prior_joint_log_likelihood(dp->hdp, fctr);
    }
    else {
        double parent_gamma = *(parent_dp->gamma);
        double log_likelihood = MINUS_INF;
        int64_t num_parent_fctrs = stSet_size(parent_dp->factors);
        Factor** parent_fctrs = (Factor**) malloc(sizeof(Factor*) * num_parent_fctrs);
        
        double next_height_unobs_log_likelihood;
        #pragma omp parallel shared(log_likelihood,next_height_unobs_log_likelihood,parent_dp,num_parent_fctrs,parent_fctrs)
        {
            #pragma omp single nowait
            next_height_unobs_log_likelihood = unobserved_factor_joint_log_likelihood(fctr, parent_dp);
            
            #pragma omp single
            {
                stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
                for (int64_t i = 0; i < num_parent_fctrs; i++) {
                    parent_fctrs[i] = (Factor*) stSet_getNext(parent_fctr_iter);
                }
                stSet_destructIterator(parent_fctr_iter);
            }
            
            double local_log_likelihood = MINUS_INF;
            double log_fctr_size;
            Factor* parent_fctr;
            
            #pragma omp for nowait
            for (int64_t i = 0; i < num_parent_fctrs; i++) {
                parent_fctr = parent_fctrs[i];
                log_fctr_size = log(stSet_size(get_factor_children(parent_fctr)));
                local_log_likelihood = add_logs(local_log_likelihood,
                                                log_fctr_size + factor_parent_joint_log_likelihood(fctr, parent_fctr));
            }
            
            #pragma omp critical
            log_likelihood = add_logs(log_likelihood, local_log_likelihood);
        }
        free(parent_fctrs);
        
        log_likelihood = add_logs(log_likelihood,
                                  log(parent_gamma) + next_height_unobs_log_likelihood);
        
        log_likelihood -= log(parent_gamma + parent_dp->num_factor_children);
        
        return log_likelihood;
    }
}

//double unobserved_factor_joint_log_likelihood(Factor* fctr, DirichletProcess* dp) {
//    DirichletProcess* parent_dp = dp->parent;
//    if (parent_dp == NULL) {
//        return prior_joint_log_likelihood(dp->hdp, fctr);
//    }
//    else {
//        double parent_gamma = *(parent_dp->gamma);
//        
//        double log_likelihood = MINUS_INF;
//        stSetIterator* parent_fctr_iter = stSet_getIterator(parent_dp->factors);
//        Factor* parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
//        double log_fctr_size;
//        while (parent_fctr != NULL) {
//            log_fctr_size = log((double) stSet_size(parent_fctr->children));
//            log_likelihood = add_logs(log_likelihood,
//                                      log_fctr_size + factor_parent_joint_log_likelihood(fctr, parent_fctr));
//            parent_fctr = (Factor*) stSet_getNext(parent_fctr_iter);
//        }
//        stSet_destructIterator(parent_fctr_iter);
//        
//        log_likelihood = add_logs(log_likelihood,
//                                  log(parent_gamma) + unobserved_factor_joint_log_likelihood(fctr, parent_dp));
//        
//        log_likelihood -= log(parent_gamma + (double) parent_dp->num_factor_children);
//        
//        return log_likelihood;
//    }
//}

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
    if (fctr->factor_type == BASE || fctr->factor_type == MIDDLE) {
        stSetIterator* child_fctr_iter = stSet_getIterator(get_factor_children(fctr));
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
HierarchicalDirichletProcess* new_hier_dir_proc(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop, int64_t sampling_grid_length, double mu,
                                                double nu, double alpha, double beta) {
    if (nu <= 0.0) {
        fprintf(stderr, "Nu parameter of Normal-Inverse-Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (alpha <= 0.0) {
        fprintf(stderr, "Alpha parameter of Normal-Inverse-Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (beta <= 0.0) {
        fprintf(stderr, "Beta parameter of Normal-Inverse-Gamma distribution must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (gamma != NULL) {
        for (int64_t i = 0; i < depth; i++) {
            if (gamma[i] <= 0) {
                fprintf(stderr, "All concentration parameters gamma must be postive.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if (num_dps < 2) {
        fprintf(stderr, "Hierarchical Dirichlet process formalism requires at least two Dirichlet processes.\n");
        exit(EXIT_FAILURE);
    }
    
    double* grid = linspace(sampling_grid_start, sampling_grid_stop, sampling_grid_length);

    HierarchicalDirichletProcess* hdp = (HierarchicalDirichletProcess*) malloc(sizeof(HierarchicalDirichletProcess));

    // normal-inverse gamma parameters

    hdp->mu = mu;
    hdp->nu = nu;
    hdp->alpha = alpha;
    hdp->beta = beta;

    hdp->gamma = gamma;
    hdp->depth = depth;

    hdp->finalized = false;
    hdp->fuzzy_assignments = false;
    hdp->has_data = false;

    hdp->num_dps = num_dps;
    DirichletProcess** dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
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
    
    hdp->data_pt_fuzzy_dp_ids = NULL;
    hdp->data_pt_fuzzy_dp_probs = NULL;
    hdp->data_pt_num_fuzzy_dps = NULL;


    hdp->sample_gamma = false;
    hdp->gamma_alpha = NULL;
    hdp->gamma_beta = NULL;
    hdp->s_aux_vector = NULL;
    hdp->w_aux_vector = NULL;
    
    hdp->distr_metric_memos = stSet_construct2(&destroy_distr_metric_memo);

    return hdp;
}

// Gamma prior on concentration parameters
HierarchicalDirichletProcess* new_hier_dir_proc_2(int64_t num_dps, int64_t depth, double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start, double sampling_grid_stop,
                                                  int64_t sampling_grid_length, double mu, double nu, double alpha,
                                                  double beta) {

    for (int64_t i = 0; i < depth; i++) {
        if (gamma_alpha[i] <= 0.0) {
            fprintf(stderr, "Alpha parameter of Gamma distribution must be positive.\n");
            exit(EXIT_FAILURE);
        }
        if (gamma_beta[i] <= 0.0) {
            fprintf(stderr, "Beta parameter of Gamma distribution must be positive.\n");
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

    for (int64_t i = 0; i < num_dps; i++) {
        w[i] = 1.0;
        s[i] = false;
    }
    
    // init to prior expected value
    double* gamma = (double*) malloc(sizeof(double) * depth);
    hdp->gamma = gamma;
    for (int64_t i = 0; i < depth; i++) {
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
    free(hdp->gamma_alpha);
    free(hdp->gamma_beta);
    free(hdp->w_aux_vector);
    free(hdp->s_aux_vector);
    stSet_destruct(hdp->distr_metric_memos);
    free(hdp);
}

void establish_base_dp(HierarchicalDirichletProcess* hdp) {
    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;

    DirichletProcess* dp;
    for (int64_t i = 0; i < num_dps; i++) {
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
    int64_t num_dps = hdp->num_dps;
    bool* visited = (bool*) malloc(sizeof(bool) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
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

void verify_tree_depth(HierarchicalDirichletProcess* hdp, DirichletProcess* dp, int64_t current_depth,
                       int64_t leaf_depth) {
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

void verify_valid_fuzzy_dp_assignments(HierarchicalDirichletProcess* hdp) {
    int64_t length = hdp->data_length;
    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    int64_t** fuzzy_dp_ids = hdp->data_pt_fuzzy_dp_ids;
    int64_t* num_fuzzy_dps = hdp->data_pt_num_fuzzy_dps;
    
    int64_t id;
    int64_t* ids;
    int64_t num_ids;
    DirichletProcess* dp;
    for (int64_t i = 0; i < length; i++) {
        ids = fuzzy_dp_ids[i];
        num_ids = num_fuzzy_dps[i];
        for (int64_t j = 0; j < num_ids; j++) {
            id = ids[j];
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
}

void verify_valid_definite_dp_assignments(HierarchicalDirichletProcess* hdp) {
    int64_t length = hdp->data_length;
    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    int64_t* dp_ids = hdp->data_pt_dp_id;

    int64_t id;
    DirichletProcess* dp;
    for (int64_t i = 0; i < length; i++) {
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

void verify_valid_dp_assignments(HierarchicalDirichletProcess* hdp) {
    if (hdp->fuzzy_assignments) {
        verify_valid_fuzzy_dp_assignments(hdp);
    }
    else {
        verify_valid_definite_dp_assignments(hdp);
    }
}

void mark_fuzzy_observed_dps(HierarchicalDirichletProcess* hdp) {
    int64_t length = hdp->data_length;
    DirichletProcess** dps = hdp->dps;
    int64_t** fuzzy_dp_ids = hdp->data_pt_fuzzy_dp_ids;
    int64_t* num_fuzzy_dps = hdp->data_pt_num_fuzzy_dps;
    int64_t grid_length = hdp->grid_length;
    
    DirichletProcess* dp;
    double* pdf;
    int64_t* ids;
    int64_t num_ids;
    for (int64_t i = 0; i < length; i++) {
        ids = fuzzy_dp_ids[i];
        num_ids = num_fuzzy_dps[i];
        for (int64_t j = 0; j < num_ids; j++) {
            dp = dps[ids[j]];
            while (dp != NULL) {
                if (dp->observed) {
                    break;
                }
                dp->observed = true;
                
                pdf = (double*) malloc(sizeof(double) * grid_length);
                dp->posterior_predictive = pdf;
                for (int64_t k = 0; k < grid_length; k++) {
                    pdf[k] = 0.0;
                }
                
                dp = dp->parent;
            }
        }
    }
}

void mark_definite_observed_dps(HierarchicalDirichletProcess* hdp) {
    int64_t length = hdp->data_length;
    DirichletProcess** dps = hdp->dps;
    int64_t* dp_ids = hdp->data_pt_dp_id;
    int64_t grid_length = hdp->grid_length;

    DirichletProcess* dp;
    double* pdf;
    int64_t id;
    for (int64_t i = 0; i < length; i++) {
        id = dp_ids[i];
        dp = dps[id];
        while (dp != NULL) {
            if (dp->observed) {
                break;
            }
            dp->observed = true;

            pdf = (double*) malloc(sizeof(double) * grid_length);
            dp->posterior_predictive = pdf;
            for (int64_t j = 0; j < grid_length; j++) {
                pdf[j] = 0.0;
            }

            dp = dp->parent;
        }
    }
}

void mark_observed_dps(HierarchicalDirichletProcess* hdp) {
    if (hdp->fuzzy_assignments) {
        mark_fuzzy_observed_dps(hdp);
    }
    else {
        mark_definite_observed_dps(hdp);
    }
}

void init_factors_internal(DirichletProcess* dp, Factor* parent_fctr, stList** data_pt_fctr_lists) {
    if (!dp->observed) {
        return;
    }
    Factor* middle_fctr = new_middle_factor(dp);
    Factor** middle_fctr_parent_ptr = get_factor_parent_ptr(middle_fctr);
    *middle_fctr_parent_ptr = parent_fctr;
    stSet_insert(get_factor_children(parent_fctr), (void*) middle_fctr);
    
    if (stList_length(dp->children) == 0) {
        stSet* middle_fctr_children = get_factor_children(middle_fctr);
        
        stListIterator* data_pt_fctr_iter = stList_getIterator(data_pt_fctr_lists[dp->id]);
        Factor* data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        while (data_pt_fctr != NULL) {
            Factor** data_pt_fctr_parent_ptr = get_factor_parent_ptr(data_pt_fctr);
            *data_pt_fctr_parent_ptr = middle_fctr;
            stSet_insert(middle_fctr_children, (void*) data_pt_fctr);
            data_pt_fctr = (Factor*) stList_getNext(data_pt_fctr_iter);
        }
        stList_destructIterator(data_pt_fctr_iter);
    }
    else {
        stListIterator* child_dp_iter = stList_getIterator(dp->children);
        DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        while (child_dp != NULL) {
            init_factors_internal(child_dp, middle_fctr, data_pt_fctr_lists);
            child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
        }
        stList_destructIterator(child_dp_iter);
    }
}

stList** get_definite_data_pt_factor_lists(HierarchicalDirichletProcess* hdp) {
    int64_t data_length = hdp->data_length;
    int64_t* data_pt_dp_id = hdp->data_pt_dp_id;
    double* data = hdp->data;
    int64_t num_dps = hdp->num_dps;
    
    stList** data_pt_fctr_lists = (stList**) malloc(sizeof(stList*) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
        data_pt_fctr_lists[i] = NULL;
    }
    
    Factor* data_pt_fctr;
    int64_t dp_id;
    stList* fctr_list;
    for (int64_t data_pt_idx = 0; data_pt_idx < data_length; data_pt_idx++) {
        dp_id = data_pt_dp_id[data_pt_idx];
        fctr_list = data_pt_fctr_lists[dp_id];
        if (fctr_list == NULL) {
            fctr_list = stList_construct();
            data_pt_fctr_lists[dp_id] = fctr_list;
        }
        data_pt_fctr = new_data_pt_factor(data[data_pt_idx]);
        stList_append(fctr_list, (void*) data_pt_fctr);
    }
    
    return data_pt_fctr_lists;
}

stList** get_fuzzy_data_pt_factor_lists(HierarchicalDirichletProcess* hdp) {
    int64_t data_length = hdp->data_length;
    int64_t* data_pt_num_fuzzy_dps = hdp->data_pt_num_fuzzy_dps;
    int64_t** data_pt_fuzzy_dp_ids = hdp->data_pt_fuzzy_dp_ids;
    double** data_pt_fuzzy_dp_probs = hdp->data_pt_fuzzy_dp_probs;
    double* data = hdp->data;
    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;
    
    stList** data_pt_fctr_lists = (stList**) malloc(sizeof(stList*) * num_dps);
    for (int64_t i = 0; i < num_dps; i++) {
        data_pt_fctr_lists[i] = NULL;
    }
    
    Factor* data_pt_fctr;
    int64_t dp_id;
    int64_t num_fuzzy_dps;
    int64_t* fuzzy_dp_ids;
    DirichletProcess* dp_assignment;
    stList* fctr_list;
    DirichletProcess** fuzzy_dps;
    double* fuzzy_dp_probs;
    
    for (int64_t data_pt_idx = 0; data_pt_idx < data_length; data_pt_idx++) {
        fuzzy_dp_ids = data_pt_fuzzy_dp_ids[data_pt_idx];
        num_fuzzy_dps = data_pt_num_fuzzy_dps[data_pt_idx];
        fuzzy_dp_probs = data_pt_fuzzy_dp_probs[data_pt_idx];
        
        fuzzy_dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_fuzzy_dps);
        for (int64_t j = 0; j < num_fuzzy_dps; j++) {
            fuzzy_dps[j] = dps[fuzzy_dp_ids[j]];
        }
        data_pt_fctr = new_fuzzy_data_pt_factor(data[data_pt_idx], fuzzy_dps, fuzzy_dp_probs, num_fuzzy_dps);
        dp_assignment = sample_fuzzy_data_pt_dp_assignment(data_pt_fctr);
        
        dp_id = dp_assignment->id;
        fctr_list = data_pt_fctr_lists[dp_id];
        if (fctr_list == NULL) {
            fctr_list = stList_construct();
            data_pt_fctr_lists[dp_id] = fctr_list;
        }
        stList_append(fctr_list, (void*) data_pt_fctr);
    }
    
    return data_pt_fctr_lists;
}

void init_factors(HierarchicalDirichletProcess* hdp) {
    
    // decentralize data into factors
    stList** data_pt_fctr_lists;
    if (hdp->fuzzy_assignments) {
        data_pt_fctr_lists = get_fuzzy_data_pt_factor_lists(hdp);
        
        free(hdp->data_pt_fuzzy_dp_probs);
        hdp->data_pt_fuzzy_dp_probs = NULL;
        
        free(hdp->data_pt_fuzzy_dp_ids);
        hdp->data_pt_fuzzy_dp_ids = NULL;

        free(hdp->data_pt_num_fuzzy_dps);
        hdp->data_pt_num_fuzzy_dps = NULL;
    }
    else {
        data_pt_fctr_lists = get_definite_data_pt_factor_lists(hdp);
        
        free(hdp->data_pt_dp_id);
        hdp->data_pt_dp_id = NULL;
    }
    free(hdp->data);
    hdp->data = NULL;

    DirichletProcess* base_dp = hdp->base_dp;
    Factor* root_factor = new_base_factor(hdp);
    
    stListIterator* child_dp_iter = stList_getIterator(base_dp->children);
    DirichletProcess* child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    while (child_dp != NULL) {
        init_factors_internal(child_dp, root_factor, data_pt_fctr_lists);
        child_dp = (DirichletProcess*) stList_getNext(child_dp_iter);
    }
    stList_destructIterator(child_dp_iter);
    
    int64_t num_dps = hdp->num_dps;
    for (int64_t i = 0; i < num_dps; i++) {
        if (data_pt_fctr_lists[i] != NULL) {
            stList_destruct(data_pt_fctr_lists[i]);
        }
    }
    free(data_pt_fctr_lists);
    
    double mean, sum_sq_devs;
    int64_t num_data;
    get_factor_stats(root_factor, &mean, &sum_sq_devs, &num_data);
    add_update_base_factor_params(root_factor, mean, sum_sq_devs, (double) num_data);
    
    int64_t fctr_child_count;
    DirichletProcess* dp;
    Factor* fctr;
    stSetIterator* fctr_iter;
    DirichletProcess** dps = hdp->dps;
    for (int64_t i = 0; i < num_dps; i++) {
        dp = dps[i];
        
        fctr_child_count = 0;

        fctr_iter = stSet_getIterator(dp->factors);
        fctr = (Factor*) stSet_getNext(fctr_iter);
        while (fctr != NULL) {
            fctr_child_count += stSet_size(get_factor_children(fctr));
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

void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int64_t child_id, int64_t parent_id) {
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

void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int64_t* dp_ids, int64_t length) {
    if (hdp->has_data) {
        fprintf(stderr, "Hierarchical Dirichlet process must be reset before passing new data.\n");
        exit(EXIT_FAILURE);
    }

    hdp->data = data;
    hdp->data_pt_dp_id = dp_ids;
    hdp->data_length = length;
    
    hdp->fuzzy_assignments = false;
    hdp->has_data = true;

    if (hdp->finalized) {
        finalize_data(hdp);
    }
}

void pass_fuzzy_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int64_t** fuzzy_dp_ids,
                            double** assignment_probs, int64_t* num_fuzzy_dps, int64_t length) {
    if (hdp->has_data) {
        fprintf(stderr, "Hierarchical Dirichlet process must be reset before passing new data.\n");
        exit(EXIT_FAILURE);
    }
    
    hdp->data = data;
    hdp->data_pt_fuzzy_dp_probs = assignment_probs;
    hdp->data_pt_fuzzy_dp_ids = fuzzy_dp_ids;
    hdp->data_pt_num_fuzzy_dps = num_fuzzy_dps;
    hdp->data_length = length;
    
    hdp->fuzzy_assignments = true;
    hdp->has_data = true;
    
    if (hdp->finalized) {
        finalize_data(hdp);
    }
}

void finalize_hdp_structure(HierarchicalDirichletProcess* hdp) {
    establish_base_dp(hdp);
    verify_dp_tree(hdp);
    verify_tree_depth(hdp, hdp->base_dp, 0, hdp->depth - 1);
    
    hdp->finalized = true;
    
    if (hdp->has_data) {
        finalize_data(hdp);
    }
}

void reset_distr_metric_memo(DistributionMetricMemo* memo) {
    int64_t num_distrs = memo->num_distrs;
    int64_t num_entries = ((num_distrs - 1) * num_distrs) / 2;
    double* memo_entries = memo->memo_matrix;
    for (int64_t i = 0; i < num_entries; i++) {
        memo_entries[i] = -1.0;
    }
}

void reset_hdp_data(HierarchicalDirichletProcess* hdp) {
    if (!hdp->has_data) {
        return;
    }
    
    hdp->data_length = 0;

    DirichletProcess** dps = hdp->dps;
    int64_t num_dps = hdp->num_dps;
    
    destroy_dir_proc_factor_tree(hdp->base_dp);
    
    DirichletProcess* dp;
    for (int64_t i = 0; i < num_dps; i++) {
        dp = dps[i];

        free(dp->posterior_predictive);
        dp->posterior_predictive = NULL;

        free(dp->spline_slopes);
        dp->spline_slopes = NULL;

        dp->observed = false;
    }
    
    stSetIterator* memo_iter = stSet_getIterator(hdp->distr_metric_memos);
    DistributionMetricMemo* memo = stSet_getNext(memo_iter);
    while (memo != NULL) {
        reset_distr_metric_memo(memo);
        memo = stSet_getNext(memo_iter);
    }
    stSet_destructIterator(memo_iter);
    
    hdp->splines_finalized = false;
    
    hdp->samples_taken = 0;
    
    if (hdp->sample_gamma) {
        double* gamma = hdp->gamma;
        double* gamma_alpha = hdp->gamma_alpha;
        double* gamma_beta = hdp->gamma_beta;
        
        for (int64_t depth = 0; depth < hdp->depth; depth++) {
            gamma[depth] = gamma_alpha[depth] / gamma_beta[depth];
        }
        
        double* w = hdp->w_aux_vector;
        bool* s = hdp->s_aux_vector;
        
        for (int64_t i = 0; i < num_dps; i++) {
            w[i] = 1.0;
            s[i] = false;
        }
    }
    
    hdp->has_data = false;
}

void unassign_from_parent(Factor* fctr) {
    //printf("getting parent factor\n");
    Factor* parent;
    switch (fctr->factor_type) {
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            fctr_data->parent = NULL;
            break;
        }
        case DATA_PT:
        {
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            fctr_data->parent = NULL;
            break;
        }
        case FUZZY_DATA_PT:
        {
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            parent = fctr_data->parent;
            fctr_data->parent = NULL;
            break;
        }
        default:
        {
            fprintf(stderr, "Cannot unassign base factor's parent.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    //printf("getting parent dp\n");
    DirichletProcess* parent_dp = get_factor_dir_proc(parent);
    
    //printf("getting base factor\n");
    Factor* base_fctr = get_base_factor(parent);
    //printf("getting base dp\n");
    DirichletProcess* base_dp = get_factor_dir_proc(base_fctr);
    
    stSet* parents_children = get_factor_children(parent);
    
    //printf("removing from set\n");
    stSet_remove(parents_children, (void*) fctr);
    //printf("decrementing count\n");
    (parent_dp->num_factor_children)--;
    
    if (stSet_size(parents_children) == 0) {
        //printf("destroying parent factor\n");
        destroy_factor(parent);
    }

    int64_t num_reassign;
    double mean_reassign;
    double sum_sq_devs;
    
    //printf("getting stats\n");
    get_factor_stats(fctr, &mean_reassign, &sum_sq_devs, &num_reassign);
    
    // check to see if base factor has been destroyed
    if (stSet_search(base_dp->factors, (void*) base_fctr) != NULL) {
        //printf("updating params\n");
        remove_update_base_factor_params(base_fctr, mean_reassign, sum_sq_devs, (double) num_reassign);
    }
    
    //printf("caching stats\n");
    if (fctr->factor_type == MIDDLE) {
        DirichletProcess* dp = get_factor_dir_proc(fctr);
        dp->cached_factor_mean = mean_reassign;
        dp->cached_factor_size = num_reassign;
        dp->cached_factor_sum_sq_dev = sum_sq_devs;
    }
}

void assign_to_parent(Factor* fctr, Factor* parent, bool update_params) {
    
    switch (fctr->factor_type) {
        case MIDDLE:
        {
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            fctr_data->parent = parent;
            break;
        }
        case DATA_PT:
        {
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            fctr_data->parent = parent;
            break;
        }
        case FUZZY_DATA_PT:
        {
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            fctr_data->parent = parent;
            break;
        }
        default:
        {
            fprintf(stderr, "Cannot assign base factor's parent.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    stSet* parents_children = get_factor_children(parent);
    stSet_insert(parents_children, (void*) fctr);
    
    DirichletProcess* parent_dp = get_factor_dir_proc(parent);
    (parent_dp->num_factor_children)++;

    if (!update_params) {
        return;
    }
    
    Factor* base_fctr = get_base_factor(parent);

    if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
        double data_pt = get_factor_data_pt(fctr);
        add_update_base_factor_params(base_fctr, data_pt, 0.0, 1.0);
    }
    else {
        DirichletProcess* dp = get_factor_dir_proc(fctr);
        add_update_base_factor_params(base_fctr, dp->cached_factor_mean, dp->cached_factor_sum_sq_dev,
                                      (double) dp->cached_factor_size);
    }
}

//Factor* sample_from_data_pt_factor(Factor* fctr, DirichletProcess* dp) {
//    if (fctr->factor_type != DATA_PT) {
//        fprintf(stderr, "Attempted a data point factor sample from non-data point factor.\n");
//        exit(EXIT_FAILURE);
//    }
//    
//    stSet* pool = dp->factors;
//    int64_t num_fctrs = stSet_size(pool);
//    
//    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);
//    
//    double* cdf = (double*) malloc(sizeof(double) * (num_fctrs + 1));
//    double cumul = 0.0;
//    
//    stSetIterator* pool_iter = stSet_getIterator(pool);
//    Factor* fctr_option;
//    double fctr_size;
//    for (int64_t i = 0; i < num_fctrs; i++) {
//        fctr_option = (Factor*) stSet_getNext(pool_iter);
//        fctr_order[i] = fctr_option;
//        
//        fctr_size = (double) stSet_size(fctr_option->children);
//        cumul += fctr_size * data_pt_factor_parent_likelihood(fctr, fctr_option);
//        cdf[i] = cumul;
//    }
//    stSet_destructIterator(pool_iter);
//    
//    double gamma_param = *(dp->gamma);
//    cumul += gamma_param * unobserved_factor_likelihood(fctr, dp);
//    cdf[num_fctrs] = cumul;
//    
//    int64_t choice_idx = bisect_left(rand_uniform(cumul), cdf, num_fctrs + 1);
//    
//    Factor* fctr_choice;
//    if (choice_idx == num_fctrs) {
//        free(fctr_order);
//        DirichletProcess* parent_dp = dp->parent;
//        if (parent_dp == NULL) {
//            fctr_choice = new_base_factor(dp->hdp);
//        }
//        else {
//            fctr_choice = new_middle_factor(dp);
//            Factor* new_fctr_parent = sample_from_data_pt_factor(fctr, parent_dp);
//            assign_to_parent(fctr_choice, new_fctr_parent, false);
//        }
//    }
//    else {
//        fctr_choice = fctr_order[choice_idx];
//        free(fctr_order);
//    }
//    
//    return fctr_choice;
//}

Factor* sample_from_data_pt_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != DATA_PT && fctr->factor_type != FUZZY_DATA_PT) {
        fprintf(stderr, "Attempted a data point factor sample from non-data point factor.\n");
        exit(EXIT_FAILURE);
    }
    
    stSet* pool = dp->factors;
    int64_t num_fctrs = stSet_size(pool);
    
    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);


    stSetIterator* pool_iter = stSet_getIterator(pool);
    for (int64_t i = 0; i < num_fctrs; i++) {
        Factor* fctr_option = (Factor*) stSet_getNext(pool_iter);
        fctr_order[i] = fctr_option;
    }
    stSet_destructIterator(pool_iter);
    
    double* probs = (double*) malloc(sizeof(double) * num_fctrs);
    double new_fctr_prob;
    #pragma omp parallel shared(new_fctr_prob,probs)
    {
        #pragma omp single nowait 
        new_fctr_prob = (*(dp->gamma)) * unobserved_factor_likelihood(fctr, dp);
        
        Factor* fctr_option;
        #pragma omp for
        for (int64_t i = 0; i < num_fctrs; i++) {
            fctr_option = fctr_order[i];
            probs[i] = stSet_size(get_factor_children(fctr_option))
                       * data_pt_factor_parent_likelihood(fctr, fctr_option);
        }
    }
    
    double* cdf = (double*) malloc(sizeof(double) * (num_fctrs + 1));
    parallel_cdf(cdf, probs, num_fctrs, 10);
    
    cdf[num_fctrs] = cdf[num_fctrs - 1] + new_fctr_prob;
    
    int64_t choice_idx = bisect_left(rand_uniform(cdf[num_fctrs]), cdf, num_fctrs + 1);
    
    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        free(fctr_order);
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_data_pt_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent, false);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
        free(fctr_order);
    }

    return fctr_choice;
}

//Factor* sample_from_middle_factor(Factor* fctr, DirichletProcess* dp) {
//    if (fctr->factor_type != MIDDLE) {
//        fprintf(stderr, "Attempted a middle factor sample from non-middle factor.\n");
//        exit(EXIT_FAILURE);
//    }
//    
//    stSet* pool = dp->factors;
//    int64_t num_fctrs = stSet_size(pool);
//    int64_t num_choices = num_fctrs + 1;
//    
//    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);
//    double* log_probs = (double*) malloc(sizeof(double) * num_choices);
//    
//    stSetIterator* pool_iter = stSet_getIterator(pool);
//    for (int64_t i = 0; i < num_fctrs; i++) {
//        Factor* fctr_option = (Factor*) stSet_getNext(pool_iter);
//        fctr_order[i] = fctr_option;
//        log_probs[i] = log((double) stSet_size(fctr_option->children))
//                       + factor_parent_joint_log_likelihood(fctr, fctr_option);
//    }
//    stSet_destructIterator(pool_iter);
//    
//    log_probs[num_fctrs] = log(*(dp->gamma)) + unobserved_factor_joint_log_likelihood(fctr, dp);
//    
//    double* cdf = (double*) malloc(sizeof(double) * num_choices);
//    double cumul = 0.0;
//    double normalizing_const = max(log_probs, num_choices);
//    
//    for (int64_t i = 0; i < num_choices; i++) {
//        cumul += exp(log_probs[i] - normalizing_const);
//        cdf[i] = cumul;
//    }
//    
//    free(log_probs);
//    
//    int64_t choice_idx = bisect_left(rand_uniform(cumul), cdf, num_choices);
//    free(cdf);
//    
//    Factor* fctr_choice;
//    if (choice_idx == num_fctrs) {
//        free(fctr_order);
//        DirichletProcess* parent_dp = dp->parent;
//        if (parent_dp == NULL) {
//            fctr_choice = new_base_factor(dp->hdp);
//        }
//        else {
//            fctr_choice = new_middle_factor(dp);
//            Factor* new_fctr_parent = sample_from_middle_factor(fctr, parent_dp);
//            assign_to_parent(fctr_choice, new_fctr_parent, false);
//        }
//    }
//    else {
//        fctr_choice = fctr_order[choice_idx];
//        free(fctr_order);
//    }
//    
//    return fctr_choice;
//}

Factor* sample_from_middle_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type != MIDDLE) {
        fprintf(stderr, "Attempted a middle factor sample from non-middle factor.\n");
        exit(EXIT_FAILURE);
    }
    
    stSet* pool = dp->factors;
    int64_t num_fctrs = stSet_size(pool);
    int64_t num_choices = num_fctrs + 1;

    Factor** fctr_order = (Factor**) malloc(sizeof(Factor*) * num_fctrs);
    double* log_probs = (double*) malloc(sizeof(double) * num_choices);
    
    stSetIterator* pool_iter = stSet_getIterator(pool);
    for (int64_t i = 0; i < num_fctrs; i++) {
        fctr_order[i] = (Factor*) stSet_getNext(pool_iter);
    }
    stSet_destructIterator(pool_iter);
    
    
    double log_gamma_param = log(*(dp->gamma));
    double new_fctr_log_prob;
    #pragma omp parallel shared(new_fctr_log_prob,log_probs,log_gamma_param)
    {
        #pragma omp single nowait
        new_fctr_log_prob = log_gamma_param + unobserved_factor_joint_log_likelihood(fctr, dp);
        
        #pragma omp for
        for (int64_t i = 0; i < num_fctrs; i++) {
            Factor* fctr_option = fctr_order[i];
            log_probs[i] = log((double) stSet_size(get_factor_children(fctr_option)))
                           + factor_parent_joint_log_likelihood(fctr, fctr_option);
        }
    }

    log_probs[num_fctrs] = new_fctr_log_prob;
    
    double normalizing_const = parallel_max(log_probs, num_choices);
    
    parallel_add(-normalizing_const, log_probs, num_choices);
    parallel_exp(log_probs, num_choices);
    
    double* cdf = (double*) malloc(sizeof(double) * num_choices);
    parallel_cdf(cdf, log_probs, num_choices, 10);
    free(log_probs);

    int64_t choice_idx = bisect_left(rand_uniform(cdf[num_fctrs]), cdf, num_choices);
    free(cdf);

    Factor* fctr_choice;
    if (choice_idx == num_fctrs) {
        free(fctr_order);
        DirichletProcess* parent_dp = dp->parent;
        if (parent_dp == NULL) {
            fctr_choice = new_base_factor(dp->hdp);
        }
        else {
            fctr_choice = new_middle_factor(dp);
            Factor* new_fctr_parent = sample_from_middle_factor(fctr, parent_dp);
            assign_to_parent(fctr_choice, new_fctr_parent, false);
        }
    }
    else {
        fctr_choice = fctr_order[choice_idx];
        free(fctr_order);
    }

    return fctr_choice;
}

Factor* sample_factor(Factor* fctr, DirichletProcess* dp) {
    if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
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
    //printf("getting parent dp\n");
    DirichletProcess* parent_dp;
    if (fctr->factor_type != FUZZY_DATA_PT) {
        parent_dp = get_factor_dir_proc(get_factor_parent(fctr));
    }
    else {
        parent_dp = sample_fuzzy_data_pt_dp_assignment(fctr);
    }
    
    //printf("unassigning from parent\n");
    unassign_from_parent(fctr);
    //printf("sampling new parent\n");
    Factor* new_parent = sample_factor(fctr, parent_dp);
    //printf("assigning to new parent\n");
    assign_to_parent(fctr, new_parent, true);
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
    DirichletProcess* dp = get_factor_dir_proc(fctr);
    
    double gamma_param = *(dp->gamma);
    double total_children = (double) dp->num_factor_children;
    double num_fctr_children = (double) stSet_size(get_factor_children(fctr));
    double wt = num_fctr_children / (gamma_param + total_children);
    dp->base_factor_wt += wt;
    
    if (stList_length(dp->children) > 0) {
        stSetIterator* child_fctr_iter = stSet_getIterator(get_factor_children(fctr));
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

void push_factor_distr(DirichletProcess* dp, double* distr, int64_t length) {
    double* sample_collector = dp->posterior_predictive;
    double wt = dp->base_factor_wt;

    for (int64_t i = 0; i < length; i++) {
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
    int64_t length = hdp->grid_length;
    double* pdf = (double*) malloc(sizeof(double) * length);
    
    //printf("number of base factors: %"PRId64", num data: %"PRId64"\n", stSet_size(base_dp->factors), hdp->data_length);
    stSetIterator* base_fctr_iter = stSet_getIterator(base_dp->factors);
    Factor* base_fctr = (Factor*) stSet_getNext(base_fctr_iter);
    while (base_fctr != NULL) {
        cache_base_factor_weight(base_fctr);
        evaluate_posterior_predictive(base_fctr, grid, pdf, length);
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
    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess** shuffled_dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_dps);
    int64_t pos;
    for (int64_t i = 0; i < num_dps; i++) {
        pos = rand() % (i + 1);
        shuffled_dps[i] = shuffled_dps[pos];
        shuffled_dps[pos] = dps[i];
    }
    return shuffled_dps;
}

void sample_dp_factors(DirichletProcess* dp, int64_t* iter_counter, int64_t burn_in, int64_t thinning,
                       int64_t* sample_counter, int64_t num_samples) {
    
    if (!dp->observed) {
        //printf("dp %"PRId64" is unobserved, not sampling\n", dp->id);
        return;
    }
    //printf("sampling factors from dp %"PRId64"\n", dp->id);
    int64_t iter = *iter_counter;
    int64_t samples_taken = *sample_counter;
    
    // have to pre-allocate the array of sampling factors in case reassignment triggers
    // destruction of the set the iterator is iterating through
    int64_t num_factor_children = dp->num_factor_children;
    Factor** sampling_fctrs = (Factor**) malloc(sizeof(Factor*) * num_factor_children);
    int64_t i = 0;
    
    //printf("allocating array of factors\n");
    stSetIterator* fctr_iter = stSet_getIterator(dp->factors);
    Factor* fctr = (Factor*) stSet_getNext(fctr_iter);
    stSetIterator* child_fctr_iter;
    Factor* child_fctr;
    while (fctr != NULL) {
        child_fctr_iter = stSet_getIterator(get_factor_children(fctr));
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
    
    //printf("sampling\n");
    for (int64_t j = 0; j < num_factor_children; j++) {
        //printf("sampling factor %"PRId64"\n", j);
        gibbs_factor_iteration(sampling_fctrs[j]);
        iter++;
        
        if (iter % thinning == 0) {
            
            if (iter > burn_in) {
                //printf("taking distribution sample\n");
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
    //printf("finished sampling from dp %"PRId64"\n", dp->id);
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
    int64_t num_dps = hdp->num_dps;
    DirichletProcess* dp;
    for (int64_t id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        w[id] = sample_auxilliary_w(dp);
        s[id] = sample_auxilliary_s(dp);
    }
}

void sample_base_gamma_internal(HierarchicalDirichletProcess* hdp, double log_w, int64_t num_factors) {
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

void sample_middle_gammas_internal(HierarchicalDirichletProcess* hdp, int64_t depth,
                         double sum_log_w, int64_t sum_s, int64_t num_depth_fctrs) {
    double gamma_alpha = hdp->gamma_alpha[depth];
    double gamma_beta = hdp->gamma_beta[depth];

    float gamma_alpha_post = (float) (gamma_alpha + (double) (num_depth_fctrs - sum_s));
    float gamma_beta_post = (float) (gamma_beta - sum_log_w);
    // note: different parameterization switches alpha and beta
    hdp->gamma[depth] = (double) gengam(gamma_beta_post, gamma_alpha_post);
}

void sample_gammas(HierarchicalDirichletProcess* hdp, int64_t* iter_counter, int64_t burn_in,
                   int64_t thinning, int64_t* sample_counter, int64_t num_samples) {
    int64_t iter = *iter_counter;
    int64_t samples_taken = *sample_counter;

    int64_t tree_depth = hdp->depth;
    double* w = hdp->w_aux_vector;
    bool* s = hdp->s_aux_vector;

    int64_t* num_depth_fctrs = (int64_t*) malloc(sizeof(int64_t) * tree_depth);
    double* sum_log_w = (double*) malloc(sizeof(double) * tree_depth);
    int64_t* sum_s = (int64_t*) malloc(sizeof(int64_t) * tree_depth);

    for (int64_t depth = 0; depth < tree_depth; depth++) {
        num_depth_fctrs[depth] = 0;
        sum_log_w[depth] = 0.0;
        sum_s[depth] = 0;
    }

    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    int64_t dp_depth;
    for (int64_t id = 0; id < num_dps; id++) {
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }
        dp_depth = dp->depth;
        num_depth_fctrs[dp_depth] += stSet_size(dp->factors);
        sum_log_w[dp_depth] += log(w[id]);
        if (s[id]) sum_s[dp_depth]++;
    }

    for (int64_t depth = 0; depth < tree_depth; depth++) {
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

void sample_gamma_params(HierarchicalDirichletProcess* hdp, int64_t* iter_counter, int64_t burn_in,
                         int64_t thinning, int64_t* sample_counter, int64_t num_samples) {
    sample_gamma_aux_vars(hdp);
    sample_gammas(hdp, iter_counter, burn_in, thinning, sample_counter, num_samples);
}

double snapshot_joint_log_density_internal(Factor* fctr) {
    if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
        return log(data_pt_factor_parent_likelihood(fctr, get_factor_parent(fctr)));
    }
    else {
        double log_density = 0.0;
        stSetIterator* child_fctr_iter = stSet_getIterator(get_factor_children(fctr));
        Factor* child_fctr = stSet_getNext(child_fctr_iter);
        while (child_fctr != NULL) {
            log_density += snapshot_joint_log_density_internal(child_fctr);
            child_fctr = stSet_getNext(child_fctr_iter);
        }
        stSet_destructIterator(child_fctr_iter);
        return log_density;
    }
}

double snapshot_joint_log_density(HierarchicalDirichletProcess* hdp) {
    double log_density = 0.0;
    stSetIterator* base_fctr_iter = stSet_getIterator(hdp->base_dp->factors);
    Factor* base_fctr = stSet_getNext(base_fctr_iter);
    while (base_fctr != NULL) {
        log_density += snapshot_joint_log_density_internal(base_fctr);
        base_fctr = stSet_getNext(base_fctr_iter);
    }
    stSet_destructIterator(base_fctr_iter);
    return log_density;
}


int64_t* snapshot_num_factors(HierarchicalDirichletProcess* hdp, int64_t* length_out) {
    int64_t length = hdp->num_dps;
    *length_out = length;
    int64_t* snapshot = (int64_t*) malloc(sizeof(int64_t) * length);
    
    DirichletProcess** dps = hdp->dps;
    for (int64_t i = 0; i < length; i++) {
        snapshot[i] = (int64_t) stSet_size((dps[i])->factors);
    }
    return snapshot;
}

double* snapshot_gamma_params(HierarchicalDirichletProcess* hdp, int64_t* length_out) {
    int64_t length = hdp->depth;
    *length_out = length;
    
    double* snapshot = (double*) malloc(sizeof(double) * length);
    double* gammas = hdp->gamma;
    for (int64_t i = 0; i < length; i++) {
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
    else if (fctr->factor_type == DATA_PT || fctr->factor_type == FUZZY_DATA_PT) {
        Factor* parent_fctr = get_factor_parent(fctr);
        DirichletProcess* parent_dp = get_factor_dir_proc(parent_fctr);
        stSet* pool = parent_dp->factors;
        int64_t num_fctrs = stSet_size(pool);
        
        stSetIterator* pool_iter = stSet_getIterator(pool);
        Factor* fctr_option;
        double fctr_size;
        double prob;
        for (int64_t i = 0; i < num_fctrs; i++) {
            fctr_option = (Factor*) stSet_getNext(pool_iter);
            fctr_size = (double) stSet_size(get_factor_children(fctr_option));
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
        DirichletProcess* dp = get_factor_dir_proc(fctr);
        DirichletProcess* parent_dp = dp->parent;
        Factor* parent_fctr = get_factor_parent(fctr);
        
        double mean, sum_sq_devs;
        int64_t num_data;
        get_factor_stats(fctr, &mean, &sum_sq_devs, &num_data);
        
        dp->cached_factor_mean = mean;
        dp->cached_factor_size = num_data;
        dp->cached_factor_sum_sq_dev = sum_sq_devs;
        
        stSet* pool = parent_dp->factors;
        int64_t num_fctrs = stSet_size(pool);
        int64_t num_choices = num_fctrs + 1;
        
        double* log_probs = (double*) malloc(sizeof(double) * num_choices);
        
        stSetIterator* pool_iter = stSet_getIterator(pool);
        Factor* fctr_option;
        double log_prob;
        double parent_log_prob;
        double fctr_size;
        for (int64_t i = 0; i < num_fctrs; i++) {
            fctr_option = (Factor*) stSet_getNext(pool_iter);
            fctr_size = (double) stSet_size(get_factor_children(fctr_option));
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
        
        for (int64_t i = 0; i < num_choices; i++) {
            cumul += exp(log_probs[i] - normalizing_const);;
        }
        
        free(log_probs);
    }
    
    // TODO: this is a hack, makes it very inaccurate for the early iterations
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
        child_fctr_iter = stSet_getIterator(get_factor_children(fctr));
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
    
    int64_t num_dps = hdp->num_dps;
    
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    for (int64_t id = 0; id < num_dps; id++){
        dp = dps[id];
        
        if (!dp->observed) {
            continue;
        }
        
        log_likelihood += snapshot_dir_proc_log_likelihood(dp);
    }
    
    return log_likelihood;
}

void take_snapshot(HierarchicalDirichletProcess* hdp, int64_t** num_dp_fctrs_out, int64_t* num_dps_out,
                   double** gamma_params_out, int64_t* num_gamma_params_out, double* log_likelihood_out,
                   double* log_density_out) {
    
    *num_dp_fctrs_out = snapshot_num_factors(hdp, num_dps_out);
    *gamma_params_out = snapshot_gamma_params(hdp, num_gamma_params_out);
    *log_likelihood_out = snapshot_log_likelihood(hdp);
    *log_density_out = snapshot_joint_log_density(hdp);
    
}

void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int64_t num_samples, int64_t burn_in,
                            int64_t thinning, bool verbose) {
    
    execute_gibbs_sampling_with_snapshots(hdp, num_samples, burn_in, thinning, NULL, NULL, verbose);
}

void execute_gibbs_sampling_with_snapshots(HierarchicalDirichletProcess* hdp, int64_t num_samples, int64_t burn_in, int64_t thinning,
                                           void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                           void* snapshot_func_args, bool verbose) {
    if (!hdp->has_data) {
        fprintf(stderr, "Cannot perform Gibbs sampling before passing data to HDP.\n");
        exit(EXIT_FAILURE);
    }
    
    if (!hdp->finalized) {
        fprintf(stderr, "Cannot perform Gibbs sampling before finalizing HDP structure.\n");
        exit(EXIT_FAILURE);
    }
    
    int64_t prev_sweep_iter_count = 0;
    int64_t sweep_counter = 1;
    int64_t iter_counter = 0;
    int64_t sample_counter = 0;
    int64_t num_dps = hdp->num_dps;
    int64_t non_data_pt_samples = 0;

    DirichletProcess** sampling_dps;
    while (sample_counter < num_samples) {
        
        if (verbose) {
            if (sweep_counter > 1) {
                non_data_pt_samples = iter_counter - prev_sweep_iter_count - hdp->data_length;
            }
            fprintf(stderr, "Beginning sweep %"PRId64". Performed %"PRId64" sampling iterations. Previous sweep sampled from ~%"PRId64" non-data point factors. Collected %"PRId64" of %"PRId64" distribution samples.\n", sweep_counter, iter_counter, non_data_pt_samples, sample_counter, num_samples);
            prev_sweep_iter_count = iter_counter;
            sweep_counter++;
        }
        
        if (snapshot_func != NULL) {
            snapshot_func(hdp, snapshot_func_args);
        }
        //printf("get shuffled dps\n");
        sampling_dps = get_shuffled_dps(hdp);
        
        for (int64_t i = 0; i < num_dps; i++) {
            //printf("sampling dp factors\n");
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
    int64_t grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;

    int64_t num_dps = hdp->num_dps;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    double* distr;
    for (int64_t id = 0; id < num_dps; id++){
        dp = dps[id];
        if (!dp->observed) {
            continue;
        }

        distr = dp->posterior_predictive;

        for (int64_t i = 0; i < grid_length; i++) {
            distr[i] = distr[i] * inv_sample_size;
        }

        dp->spline_slopes = spline_knot_slopes(grid, distr, grid_length);
    }
    
    hdp->splines_finalized = true;
}

double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int64_t dp_id) {
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

double get_dir_proc_distance(DistributionMetricMemo* memo, int64_t dp_id_1, int64_t dp_id_2) {
    int64_t num_dps = memo->num_distrs;
    if (dp_id_1 < 0 || dp_id_2 < 0 || dp_id_1 >= num_dps || dp_id_2 >= num_dps) {
        fprintf(stderr, "Invalid Dirichlet process ID.\n");
        exit(EXIT_FAILURE);
    }
    
    if (dp_id_1 == dp_id_2) {
        return 0.0;
    }
    
    if (dp_id_1 < dp_id_2) {
        return get_dir_proc_distance(memo, dp_id_2, dp_id_1);
    }
    
    int64_t idx = ((dp_id_1 - 1) * dp_id_1) / 2 + dp_id_2;
    double* matrix = memo->memo_matrix;
    if (matrix[idx] < 0) {
        matrix[idx] = memo->metric_func(memo->hdp, dp_id_1, dp_id_2);
    }
    
    return matrix[idx];
}



double dir_proc_distance(HierarchicalDirichletProcess* hdp, int64_t dp_id_1, int64_t dp_id_2,
                         double (*dist_func)(double*, double*, double*, int64_t)) {
    if (!hdp->splines_finalized) {
        fprintf(stderr, "Cannot compute a Shannon-Jensen divergence before finalizing distributions.\n");
        exit(EXIT_FAILURE);
    }
    
    int64_t grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;
    
    
    DirichletProcess* dp_1 = hdp->dps[dp_id_1];
    DirichletProcess* dp_2 = hdp->dps[dp_id_2];
    while (!dp_1->observed) {
        dp_1 = dp_1->parent;
    }
    while (!dp_2->observed) {
        dp_2 = dp_2->parent;
    }
    
    double* distr_1 = dp_1->posterior_predictive;
    double* distr_2 = dp_2->posterior_predictive;
    
    return dist_func(grid, distr_1, distr_2, grid_length);
}

double kl_divergence(double* x, double* distr_1, double* distr_2, int64_t length) {

    double divergence = 0.0;
    double left_pt = distr_1[0] * log(distr_1[0] / distr_2[0]) + distr_2[0] * log(distr_2[0] / distr_1[0]);
    double right_pt;
    double dx;
    for (int64_t i = 1; i < length; i++) {
        right_pt = distr_1[i] * log(distr_1[i] / distr_2[i]) + distr_2[i] * log(distr_2[i] / distr_1[i]);
        
        dx = x[i] - x[i - 1];
        
        divergence += 0.5 * (left_pt + right_pt) * dx;
        
        left_pt = right_pt;
    }
    
    return divergence;
}

double dir_proc_kl_divergence(HierarchicalDirichletProcess* hdp, int64_t dp_id_1, int64_t dp_id_2) {
    return dir_proc_distance(hdp, dp_id_1, dp_id_2, &kl_divergence);
}

DistributionMetricMemo* new_kl_divergence_memo(HierarchicalDirichletProcess* hdp) {
    return new_distr_metric_memo(hdp, &dir_proc_kl_divergence);
}

double hellinger_distance(double* x, double* distr_1, double* distr_2, int64_t length) {

    double integral = 0.0;
    double left_pt = sqrt(distr_1[0] * distr_2[0]);
    double right_pt;
    double dx;
    for (int64_t i = 1; i < length; i++) {
        right_pt = sqrt(distr_1[i] * distr_2[i]);
        
        dx = x[i] - x[i - 1];
        
        integral += 0.5 * (left_pt + right_pt) * dx;
        
        left_pt = right_pt;
    }
    
    //printf("integral = %lf\n", integral);
    return sqrt(1.0 - integral);
}

double dir_proc_hellinger_distance(HierarchicalDirichletProcess* hdp, int64_t dp_id_1, int64_t dp_id_2) {
    return dir_proc_distance(hdp, dp_id_1, dp_id_2, &hellinger_distance);
}

DistributionMetricMemo* new_hellinger_distance_memo(HierarchicalDirichletProcess* hdp) {
    return new_distr_metric_memo(hdp, &dir_proc_hellinger_distance);
}

double l2_distance(double* x, double* distr_1, double* distr_2, int64_t length) {
    double integral = 0.0;
    double diff = distr_1[0] - distr_2[0];
    double left_pt = diff * diff;
    double right_pt;
    double dx;
    for (int64_t i = 1; i < length; i++) {
        diff = distr_1[i] - distr_2[i];
        right_pt = diff * diff;
        
        dx = x[i] - x[i - 1];
        
        integral += 0.5 * (left_pt + right_pt) * dx;
        
        left_pt = right_pt;
    }
    
    return sqrt(integral);
}
    
double dir_proc_l2_distance(HierarchicalDirichletProcess* hdp, int64_t dp_id_1, int64_t dp_id_2) {
    return dir_proc_distance(hdp, dp_id_1, dp_id_2, &l2_distance);
}

DistributionMetricMemo* new_l2_distance_memo(HierarchicalDirichletProcess* hdp) {
    return new_distr_metric_memo(hdp, &dir_proc_l2_distance);
}

double shannon_jensen_distance(double* x, double* distr_1, double* distr_2, int64_t length) {
    double divergence = 0.0;
    
    double mean_distr_pt = 0.5 * (distr_1[0] + distr_2[0]);
    double left_pt = 0.5 * (distr_1[0] * log(distr_1[0] / mean_distr_pt) + distr_2[0] * log(distr_2[0] / mean_distr_pt));
    double right_pt;
    double dx;
    for (int64_t i = 1; i < length; i++) {
        mean_distr_pt = 0.5 * (distr_1[i] + distr_2[i]);
        right_pt = 0.5 * (distr_1[i] * log(distr_1[i] / mean_distr_pt) + distr_2[i] * log(distr_2[i] / mean_distr_pt));
        
        dx = x[i] - x[i - 1];
        
        divergence += 0.5 * (left_pt + right_pt) * dx;
        
        left_pt = right_pt;
    }
    
    return sqrt(divergence);
}

double dir_proc_shannon_jensen_distance(HierarchicalDirichletProcess* hdp, int64_t dp_id_1, int64_t dp_id_2) {
    return dir_proc_distance(hdp, dp_id_1, dp_id_2, &shannon_jensen_distance);
}

DistributionMetricMemo* new_shannon_jensen_distance_memo(HierarchicalDirichletProcess* hdp) {
    return new_distr_metric_memo(hdp, &dir_proc_shannon_jensen_distance);
}

double dir_proc_expected_val(HierarchicalDirichletProcess* hdp, int64_t dp_id) {
    double* grid = hdp->sampling_grid;
    int64_t grid_length = hdp->grid_length;
    double* distr = hdp->dps[dp_id]->posterior_predictive;
    
    double expected_val = 0.0;
    double dx;
    for (int64_t i = 1; i < grid_length; i++) {
        dx = grid[i] - grid[i-1];
        expected_val += grid[i] * distr[i] * dx;
    }
    return expected_val;
}

double dir_proc_variance(HierarchicalDirichletProcess* hdp, int64_t dp_id) {
    double* grid = hdp->sampling_grid;
    int64_t grid_length = hdp->grid_length;
    double* distr = hdp->dps[dp_id]->posterior_predictive;
    
    double expected_val = dir_proc_expected_val(hdp, dp_id);
    double variance = 0.0;
    double dev;
    double dx;
    for (int64_t i = 1; i < grid_length; i++) {
        dx = grid[i] - grid[i-1];
        dev = grid[i] - expected_val;
        variance += dev * dev * distr[i] * dx;;
    }
    return variance;
}

double compare_hdp_distrs(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1, // this HDP is the master for grid samples
                          HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2,
                          double (*dist_func)(double*, double*, double*, int64_t)) {
    
    if (!hdp_1->splines_finalized || !hdp_2->splines_finalized) {
        fprintf(stderr, "Must finalize distributions of both hierarchical Dirichlet processes before comparing.\n");
        exit(EXIT_FAILURE);
    }
    
    int64_t num_dps_1 = hdp_1->num_dps;
    int64_t num_dps_2 = hdp_2->num_dps;
    if (dp_id_1 < 0 || dp_id_2 < 0 || dp_id_1 >= num_dps_1 || dp_id_2 >= num_dps_2) {
        fprintf(stderr, "Invalid Dirchlet process ID.\n");
        exit(EXIT_FAILURE);
    }
    
    double* grid = hdp_1->sampling_grid;
    int64_t grid_length = hdp_1->grid_length;
    
    DirichletProcess* dp_1 = hdp_1->dps[dp_id_1];
    while (!dp_1->observed) {
        dp_1 = dp_1->parent;
    }
    
    double* distr_1 = dp_1->posterior_predictive;
    
    double* distr_2 = (double*) malloc(sizeof(double) * grid_length);
    
    for (int64_t i = 0; i < grid_length; i++) {
        distr_2[i] = dir_proc_density(hdp_2, grid[i], dp_id_2);
    }
    
    return dist_func(grid, distr_1, distr_2, grid_length);
}

double compare_hdp_distrs_kl_divergence(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                        HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2) {
    
    return compare_hdp_distrs(hdp_1, dp_id_1, hdp_2, dp_id_2, &kl_divergence);
}

double compare_hdp_distrs_l2_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                        HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2) {
    
    return compare_hdp_distrs(hdp_1, dp_id_1, hdp_2, dp_id_2, &l2_distance);
}

double compare_hdp_distrs_shannon_jensen_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                        HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2) {
    
    return compare_hdp_distrs(hdp_1, dp_id_1, hdp_2, dp_id_2, &shannon_jensen_distance);
}

double compare_hdp_distrs_hellinger_distance(HierarchicalDirichletProcess* hdp_1, int64_t dp_id_1,
                                        HierarchicalDirichletProcess* hdp_2, int64_t dp_id_2) {
    
    return compare_hdp_distrs(hdp_1, dp_id_1, hdp_2, dp_id_2, &hellinger_distance);
}

void serialize_factor_tree_internal(FILE* out, Factor* fctr, int64_t parent_id, int64_t* next_fctr_id, uintptr_t data_start) {
    int64_t id = *next_fctr_id;
    (*next_fctr_id)++;
    // factor type
    switch (fctr->factor_type) {
        case BASE:
            fprintf(out, "0\t");
            break;
        case MIDDLE:
            fprintf(out, "1\t");
            break;
        case DATA_PT:
            fprintf(out, "2\t");
            break;
        case FUZZY_DATA_PT:
            fprintf(out, "3\t");
            break;
        default:
            fprintf(stderr, "Unsupported factor type.\n");
            exit(EXIT_FAILURE);
    }
    // parent id
    if (fctr->factor_type == BASE) {
        fprintf(out, "-\t");
    }
    else {
        fprintf(out, "%"PRId64"\t", parent_id);
    }
    // extra data based on type
    switch (fctr->factor_type) {
        case BASE:
        {
            // posterior params
            BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
            fprintf(out, "%.17lg;%.17lg;%.17lg;%.17lg;%.17lg", fctr_data->mu, fctr_data->nu,
                    fctr_data->alpha, fctr_data->beta, fctr_data->log_posterior_term);
            break;
        }
        case MIDDLE:
        {
            // dp id
            MiddleFactorData* fctr_data = (MiddleFactorData*) fctr->factor_data;
            fprintf(out, "%"PRId64, fctr_data->dp->id);
            break;
        }
        case DATA_PT:
        {
            // data value
            DataPtFactorData* fctr_data = (DataPtFactorData*) fctr->factor_data;
            fprintf(out, "%.17lg", fctr_data->data_pt);
            break;
        }
        case FUZZY_DATA_PT:
        {
            // data value and fuzzy assignments
            FuzzyDataPtFactorData* fctr_data = (FuzzyDataPtFactorData*) fctr->factor_data;
            
            fprintf(out, "%.17lg:", fctr_data->data_pt);
            
            int64_t num_fuzzy_dps = fctr_data->num_fuzzy_dps;
            for (int64_t i = 0; i < num_fuzzy_dps - 1; i++) {
                fprintf(out, "%"PRId64";", fctr_data->fuzzy_dps[i]->id);
            }
            fprintf(out, "%"PRId64":", fctr_data->fuzzy_dps[num_fuzzy_dps - 1]->id);
            
            for (int64_t i = 0; i < num_fuzzy_dps - 1; i++) {
                fprintf(out, "%.17lg;", fctr_data->fuzzy_dp_cdf[i]);
            }
            fprintf(out, "%.17lg", fctr_data->fuzzy_dp_cdf[num_fuzzy_dps - 1]);
            break;
        }
        default:
        {
            fprintf(stderr, "Unsupported factor type for serialization.\n");
            exit(EXIT_FAILURE);
        }
    }

    fprintf(out, "\n");
    
    if (fctr->factor_type != DATA_PT && fctr->factor_type != FUZZY_DATA_PT) {
        stSetIterator* iter = stSet_getIterator(get_factor_children(fctr));
        Factor* child_fctr = (Factor*) stSet_getNext(iter);
        while (child_fctr != NULL) {
            serialize_factor_tree_internal(out, child_fctr, id, next_fctr_id, data_start);
            child_fctr = (Factor*) stSet_getNext(iter);
        }
        stSet_destructIterator(iter);
    }
}

void serialize_hdp(HierarchicalDirichletProcess* hdp, FILE* out) {
    
    int64_t num_dps = hdp->num_dps;
//    int64_t num_data = hdp->data_length;
//    double* data = hdp->data;
//    int64_t* dp_ids = hdp->data_pt_dp_id;
//    int64_t** data_pt_fuzzy_dp_ids = hdp->data_pt_fuzzy_dp_ids;
//    double** data_pt_fuzzy_dp_probs = hdp->data_pt_fuzzy_dp_probs;
//    int64_t* data_pt_num_fuzzy_dps = hdp->data_pt_num_fuzzy_dps;
    int64_t grid_length = hdp->grid_length;
    double* grid = hdp->sampling_grid;
    int64_t depth = hdp->depth;
    double*  gamma_params = hdp->gamma;
    double* gamma_alpha = hdp->gamma_alpha;
    double* gamma_beta = hdp->gamma_beta;
    double* w_aux_vector = hdp->w_aux_vector;
    bool* s_aux_vector = hdp->s_aux_vector;
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* base_dp = hdp->base_dp;
    bool has_data = hdp->has_data;
    bool fuzzy_assignments = hdp->fuzzy_assignments;
    
    if (!hdp->finalized) {
        fprintf(stderr, "Can only serialize HierarchicalDirichletProcess with finalized structure");
        exit(EXIT_FAILURE);
    }
    // splines finalized
    fprintf(out, "%"PRId64"\n", (int64_t) hdp->splines_finalized);
    // has data
    fprintf(out, "%"PRId64"\n", (int64_t) has_data);
    // fuzzy assignments
    fprintf(out, "%"PRId64"\n", (int64_t) fuzzy_assignments);
    // sample gamma
    fprintf(out, "%"PRId64"\n", (int64_t) hdp->sample_gamma);
    // num dps
    fprintf(out, "%"PRId64"\n", num_dps);
    // num data
    fprintf(out, "%"PRId64"\n", hdp->data_length);
    // don't need to get data now that it's decentralized and structure is required to be finalized
//    // data
//    if (has_data) {
//        for (int64_t i = 0; i < num_data - 1; i++) {
//            fprintf(out, "%.17lg\t", data[i]);
//        }
//        fprintf(out, "%.17lg\n", data[num_data - 1]);
//        // dp ids
//        for (int64_t i = 0; i < num_data - 1; i++) {
//            fprintf(out, "%"PRId64"\t", dp_ids[i]);
//        }
//        fprintf(out, "%"PRId64"\n", dp_ids[num_data - 1]);
//    }
    // base params
    fprintf(out, "%.17lg\t%.17lg\t%.17lg\t%.17lg\n", hdp->mu, hdp->nu, hdp->alpha, hdp->beta);
    // sampling grid
    fprintf(out, "%.17lg\t%.17lg\t%"PRId64"\n", grid[0], grid[grid_length - 1], grid_length);
    // gamma
    for (int64_t i = 0; i < depth - 1; i++) {
        fprintf(out, "%.17lg\t", gamma_params[i]);
    }
    fprintf(out, "%.17lg\n", gamma_params[depth - 1]);
    // gamma distr params
    
    if (hdp->sample_gamma) {
        // alpha
        for (int64_t i = 0; i < depth - 1; i++) {
            fprintf(out, "%.17lg\t", gamma_alpha[i]);
        }
        fprintf(out, "%.17lg\n", gamma_alpha[depth - 1]);
        // beta
        for (int64_t i = 0; i < depth - 1; i++) {
            fprintf(out, "%.17lg\t", gamma_beta[i]);
        }
        fprintf(out, "%.17lg\n", gamma_beta[depth - 1]);
        // w
        for (int64_t i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%.17lg\t", w_aux_vector[i]);
        }
        fprintf(out, "%.17lg\n", w_aux_vector[num_dps - 1]);
        // s
        for (int64_t i = 0; i < num_dps - 1; i++) {
            fprintf(out, "%"PRId64"\t", (int64_t) s_aux_vector[i]);
        }
        fprintf(out, "%"PRId64"\n", (int64_t) s_aux_vector[num_dps - 1]);
    }
    // dp parents
    DirichletProcess* dp;
    for (int64_t i = 0; i < num_dps; i++) {
        dp = dps[i];
        
        // parent
        if (dp == base_dp) {
            fprintf(out, "-\t%"PRId64"\n", dp->num_factor_children);
        }
        else {
            fprintf(out, "%"PRId64"\t%"PRId64"\n", dp->parent->id, dp->num_factor_children);
        }
    }
    // post preds
    if (has_data) {
        double* post_pred;
        for (int64_t i = 0; i < num_dps; i++) {
            dp = dps[i];
            post_pred = dp->posterior_predictive;
            if (post_pred != NULL) {
                for (int64_t j = 0; j < grid_length - 1; j++) {
                    fprintf(out, "%.17lg\t", post_pred[j]);
                }
                fprintf(out, "%.17lg", post_pred[grid_length - 1]);
            }
            fprintf(out, "\n");
        }
    }
    // spline slopes
    if (hdp->splines_finalized) {
        double* slopes;
        for (int64_t i = 0; i < num_dps; i++) {
            dp = dps[i];
            slopes = dp->spline_slopes;
            if (slopes != NULL) {
                for (int64_t i = 0; i < grid_length - 1; i++) {
                    fprintf(out, "%.17lg\t", slopes[i]);
                }
                fprintf(out, "%.17lg", slopes[grid_length - 1]);
            }
            fprintf(out, "\n");
        }
    }
    // factors
    if (has_data) {
        int64_t next_fctr_id = 0;
        uintptr_t data_start = (uintptr_t) hdp->data;
        
        stSetIterator* iter = stSet_getIterator(base_dp->factors);
        Factor* fctr = (Factor*) stSet_getNext(iter);
        while (fctr != NULL) {
            serialize_factor_tree_internal(out, fctr, -1, &next_fctr_id, data_start);
            fctr = (Factor*) stSet_getNext(iter);
        }
        stSet_destructIterator(iter);
    }
}

HierarchicalDirichletProcess* deserialize_hdp(FILE* in) {
    // splines finalized
    char* end;
    char* line = stFile_getLineFromFile(in);
    bool splines_finalized = (bool) strtol(line, &end, 10);
    free(line);
    // has data
    line = stFile_getLineFromFile(in);
    bool has_data = (bool) strtol(line, &end, 10);
    free(line);
    // fuzzy assignments
    line = stFile_getLineFromFile(in);
    bool fuzzy_assignments = (bool) strtol(line, &end, 10);
    free(line);
    // sample gamma
    line = stFile_getLineFromFile(in);
    bool sample_gamma = (bool) strtol(line, &end, 10);
    free(line);
    // num dps
    line = stFile_getLineFromFile(in);
    int64_t num_dps = (int64_t) strtol(line, &end, 10);
    free(line);
    // num data
    line = stFile_getLineFromFile(in);
    int64_t data_length = (int64_t) strtol(line, &end, 10);
    free(line);
    
    stList* tokens;
//    double* data;
//    int64_t* dp_ids;
//    int64_t data_length;
//    if (has_data) {
//        // data
//        line = stFile_getLineFromFile(in);
//        tokens = stString_split(line);
//        data_length = stList_length(tokens);
//        data = (double*) malloc(sizeof(double) * data_length);
//        for (int64_t i = 0; i < data_length; i++) {
//            sscanf(stList_get(tokens, i), "%lf", &(data[i]));
//        }
//        free(line);
//        stList_destruct(tokens);
//        // dp ids
//        line = stFile_getLineFromFile(in);
//        tokens = stString_split(line);
//        dp_ids = (int64_t*) malloc(sizeof(int64_t) * data_length);
//        for (int64_t i = 0; i < data_length; i++) {
//            sscanf((char*) stList_get(tokens, i), "%"SCNd64, &(dp_ids[i]));
//        }
//        free(line);
//        stList_destruct(tokens);
//    }
    // base params
    line = stFile_getLineFromFile(in);
    double mu, nu, alpha, beta;
    sscanf(line, "%lg\t%lg\t%lg\t%lg", &mu, &nu, &alpha, &beta);
    free(line);
    // sampling grid
    line = stFile_getLineFromFile(in);
    double grid_start, grid_stop;
    int64_t grid_length;
    sscanf(line, "%lg\t%lg\t%"SCNd64, &grid_start, &grid_stop, &grid_length);
    free(line);
    // gamma
    line = stFile_getLineFromFile(in);
    tokens = stString_split(line);
    int64_t depth = stList_length(tokens);
    double* gamma_params = (double*) malloc(sizeof(double) * depth);
    for (int64_t i = 0; i < depth; i++) {
        sscanf((char*) stList_get(tokens, i), "%lf", &(gamma_params[i]));
    }
    free(line);
    stList_destruct(tokens);
    // gamma distr params
    double* gamma_alpha;
    double* gamma_beta;
    double* w;
    bool* s;
    int64_t s_int;
    if (sample_gamma) {
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        // gamma alpha
        gamma_alpha = (double*) malloc(sizeof(double) * depth);
        for (int64_t i = 0; i < depth; i++) {
            sscanf((char*) stList_get(tokens, i), "%lf", &(gamma_alpha[i]));
        }
        free(line);
        stList_destruct(tokens);
        // gamma beta
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        gamma_beta = (double*) malloc(sizeof(double) * depth);
        for (int64_t i = 0; i < depth; i++) {
            sscanf((char*) stList_get(tokens, i), "%lf", &(gamma_beta[i]));
        }
        free(line);
        stList_destruct(tokens);
        // w
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        w = (double*) malloc(sizeof(double) * num_dps);
        for (int64_t i = 0; i < num_dps; i++) {
            sscanf((char*) stList_get(tokens, i), "%lf", &(w[i]));
        }
        free(line);
        stList_destruct(tokens);
        // s
        line = stFile_getLineFromFile(in);
        tokens = stString_split(line);
        s = (bool*) malloc(sizeof(bool) * num_dps);
        for (int64_t i = 0; i < num_dps; i++) {
            sscanf((char*) stList_get(tokens, i), "%"SCNd64, &s_int);
            s[i] = (bool) s_int;
        }
        free(line);
        stList_destruct(tokens);
    }
    // construct hdp
    HierarchicalDirichletProcess* hdp;
    if (sample_gamma) {
        hdp = new_hier_dir_proc_2(num_dps, depth, gamma_alpha, gamma_beta, grid_start,
                                  grid_stop, grid_length, mu, nu, alpha, beta);
        for (int64_t i = 0; i < depth; i++) {
            hdp->gamma[i] = gamma_params[i];
        }
        free(gamma_params);
        for (int64_t i = 0; i < num_dps; i++) {
            hdp->w_aux_vector[i] = w[i];
            hdp->s_aux_vector[i] = s[i];
        }
        free(w);
        free(s);
    }
    else {
        hdp = new_hier_dir_proc(num_dps, depth, gamma_params, grid_start, grid_stop,
                                grid_length, mu, nu, alpha, beta);
    }
    
    DirichletProcess** dps = hdp->dps;
    DirichletProcess* dp;
    
    // dp parents and num children
    int64_t parent_id;
    int64_t num_factor_children;
    for (int64_t id = 0; id < num_dps; id++) {
        line = stFile_getLineFromFile(in);
        if (line[0] != '-') {
            sscanf(line, "%"SCNd64"\t%"SCNd64, &parent_id, &num_factor_children);
            set_dir_proc_parent(hdp, id, parent_id);
            (dps[id])->num_factor_children = num_factor_children;
        }
        else {
            sscanf(line, "-\t%"SCNd64, &num_factor_children);
            (dps[id])->num_factor_children = num_factor_children;        }
        free(line);
    }
    
    finalize_hdp_structure(hdp);
    
    hdp->data_length = data_length;
    hdp->has_data = has_data;
    hdp->fuzzy_assignments = fuzzy_assignments;
    // give it data
    if (has_data) {
        // note: don't use pass_data_to_hdp because want to manually init factors
//        hdp->data = data;
//        hdp->data_pt_dp_id = dp_ids;
        
//        verify_valid_dp_assignments(hdp);
//        mark_observed_dps(hdp);
        
        // post predictives
        double* post_pred;
        for (int64_t id = 0; id < num_dps; id++) {
            dp = dps[id];
            
            line = stFile_getLineFromFile(in);
            stList* tokens = stString_split(line);
            if (stList_length(tokens) != 0) {
                free(dp->posterior_predictive);
                dp->posterior_predictive = (double*) malloc(sizeof(double) * grid_length);
                post_pred = dp->posterior_predictive;
                for (int64_t i = 0; i < grid_length; i++) {
                    sscanf((char*) stList_get(tokens, i), "%lf\n", &(post_pred[i]));
                }
                // only has a posterior predictive allocated if has been observed
                dp->observed = true;
            }
            free(line);
            stList_destruct(tokens);
        }
    }
    
    double* spline_slopes;
    if (splines_finalized) {
        hdp->splines_finalized = true;
        for (int64_t id = 0; id < num_dps; id++) {
            dp = dps[id];
            line = stFile_getLineFromFile(in);
            stList* tokens = stString_split(line);
            if (stList_length(tokens) != 0) {
                spline_slopes = (double*) malloc(sizeof(double) * grid_length);
                dp->spline_slopes = spline_slopes;
                for (int64_t i = 0; i < grid_length; i++) {
                    sscanf((char*) stList_get(tokens, i), "%lf", &(spline_slopes[i]));
                }
            }
            free(line);
            stList_destruct(tokens);
            
        }
    }
    
    if (has_data) {
        char* type_str;
        char* parent_str;
        char* dp_str;
        char* data_val_str;
        char* params_str;
        int64_t type_int;
        int64_t dp_id;
        double data_val;
        int64_t parent_idx;
        stList* params_list;
        char* fuzzy_str;
        stList* fuzzy_list;
        stList* fuzzy_aux_list;
        char* fuzzy_aux_str;
        int64_t num_fuzzy_dps;
        int64_t fuzzy_dp_id;
        char* fuzzy_dp_cdf_str;
        double* fuzzy_dp_cdf;
        char* fuzzy_dp_id_str;
        DirichletProcess** fuzzy_dps;
        Factor* fctr;
        Factor* parent_fctr;
        
        
        stList* fctr_list = stList_construct();
        line = stFile_getLineFromFile(in);
        while (line != NULL) {
            
            tokens = stString_split(line);
            type_str = (char*) stList_get(tokens, 0);
            sscanf(type_str, "%"SCNd64, &type_int);
            
            switch (type_int) {
                case 0:
                    fctr = new_base_factor(hdp);
                    params_str = (char*) stList_get(tokens, 2);
                    params_list = stString_splitByString(params_str, ";");
                    BaseFactorData* fctr_data = (BaseFactorData*) fctr->factor_data;
                    sscanf((char*) stList_get(params_list, 0), "%lf", &(fctr_data->mu));
                    sscanf((char*) stList_get(params_list, 1), "%lf", &(fctr_data->nu));
                    sscanf((char*) stList_get(params_list, 2), "%lf", &(fctr_data->alpha));
                    sscanf((char*) stList_get(params_list, 3), "%lf", &(fctr_data->beta));
                    sscanf((char*) stList_get(params_list, 4), "%lf", &(fctr_data->log_posterior_term));
                    stList_destruct(params_list);
                    break;
                    
                case 1:
                    dp_str = (char*) stList_get(tokens, 2);
                    sscanf(dp_str, "%"SCNd64, &dp_id);
                    fctr = new_middle_factor(dps[dp_id]);
                    break;
                    
                case 2:
                    data_val_str = (char*) stList_get(tokens, 2);
                    sscanf(data_val_str, "%lf", &data_val);
                    fctr = new_data_pt_factor(data_val);
                    break;
                    
                case 3:
                    fuzzy_str = (char*) stList_get(tokens, 2);
                    fuzzy_list = stString_splitByString(fuzzy_str, ":");
                    
                    data_val_str = (char*) stList_get(fuzzy_list, 0);
                    sscanf(data_val_str, "%lf", &data_val);
                    
                    fuzzy_aux_str = (char*) stList_get(fuzzy_list, 1);
                    fuzzy_aux_list = stString_splitByString(fuzzy_aux_str, ";");
                    num_fuzzy_dps = stList_length(fuzzy_aux_list);
                    fuzzy_dps = (DirichletProcess**) malloc(sizeof(DirichletProcess*) * num_fuzzy_dps);
                    for (int64_t i = 0; i < num_fuzzy_dps; i++) {
                        fuzzy_dp_id_str = stList_get(fuzzy_aux_list, i);
                        sscanf(fuzzy_dp_id_str, "%"SCNd64, &fuzzy_dp_id);
                        fuzzy_dps[i] = dps[fuzzy_dp_id];
                    }
                    stList_destruct(fuzzy_aux_list);
                    
                    // constructor automatically constructs cdf from probs
                    // this is a hack to avoid non commutative floating point operations
                    fuzzy_dp_cdf = (double*) malloc(sizeof(double) * num_fuzzy_dps);
                    for (int64_t i = 0; i < num_fuzzy_dps; i++) {
                        fuzzy_dp_cdf[i] = 0.0;
                    }
                    
                    fctr = new_fuzzy_data_pt_factor(data_val, fuzzy_dps, fuzzy_dp_cdf, num_fuzzy_dps);
                    
                    fuzzy_aux_str = (char*) stList_get(fuzzy_list, 2);
                    fuzzy_aux_list = stString_splitByString(fuzzy_aux_str, ";");
                    for (int64_t i = 0; i < num_fuzzy_dps; i++) {
                        fuzzy_dp_cdf_str = stList_get(fuzzy_aux_list, i);
                        sscanf(fuzzy_dp_cdf_str, "%lf", &(fuzzy_dp_cdf[i]));
                    }
                    stList_destruct(fuzzy_aux_list);
                    
                    stList_destruct(fuzzy_list);
                    break;
                    
                default:
                    fprintf(stderr, "Deserialization error");
                    exit(EXIT_FAILURE);
                    break;
            }
            
            stList_append(fctr_list, (void*) fctr);
            
            // set parent if applicable
            parent_str = (char*) stList_get(tokens, 1);
            if (parent_str[0] != '-') {
                sscanf(parent_str, "%"SCNd64, &parent_idx);
                parent_fctr = (Factor*) stList_get(fctr_list, parent_idx);
                
                Factor** parent_fctr_ptr = get_factor_parent_ptr(fctr);
                *parent_fctr_ptr = parent_fctr;
                stSet_insert(get_factor_children(parent_fctr), (void*) fctr);
            }
            
            free(line);
            line = stFile_getLineFromFile(in);
        }
        stList_destruct(fctr_list);
    }
    
    return hdp;
}

