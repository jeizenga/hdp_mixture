//
//  nanopore_hdp.c
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

// in 0-based index
#define ALIGNMENT_KMER_COL 9
#define ALIGNMENT_SIGNAL_COL 13
#define NUM_ALIGNMENT_COLS 14

#define MODEL_ROW_HEADER_LENGTH 1
#define MODEL_MEAN_ENTRY 0
#define MODEL_ENTRY_LENGTH 5

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <inttypes.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "nanopore_hdp.h"
#include "sonLib.h"

void normal_inverse_gamma_params_from_minION(const char* model_filepath, double* mu_out, double* nu_out,
                                             double* alpha_out, double* beta_out) {
    
    FILE* model_file = fopen(model_filepath, "r");
    
    char* line = stFile_getLineFromFile(model_file);
    stList* tokens = stString_split(line);
    
    int64_t table_length = (stList_length(tokens) - MODEL_ROW_HEADER_LENGTH) / MODEL_ENTRY_LENGTH;
    double* means = (double*) malloc(sizeof(double) * table_length);
    
    int64_t offset = MODEL_ROW_HEADER_LENGTH + MODEL_MEAN_ENTRY;
    char* mean_str;
    for (int i = 0; i < table_length; i++) {
        mean_str = (char*) stList_get(tokens, offset + i * MODEL_ENTRY_LENGTH);
        sscanf(mean_str, "%lf", &(means[i]));
    }
    
    free(line);
    stList_destruct(tokens);
    
    normal_inverse_gamma_params(means, table_length, mu_out, nu_out, alpha_out, beta_out);
    
    fclose(model_file);
}

// fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int64_t sampling_grid_length,
                                         const char* model_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(model_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc(num_dps, depth, gamma, sampling_grid_start, sampling_grid_stop,
                             sampling_grid_length, mu, nu, alpha, beta);
}

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int64_t num_dps, int64_t depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int64_t sampling_grid_length,
                                           const char* model_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(model_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc_2(num_dps, depth, gamma_alpha, gamma_beta, sampling_grid_start,
                               sampling_grid_stop, sampling_grid_length, mu, nu, alpha, beta);
}

void update_hdp_from_alignment(HierarchicalDirichletProcess* hdp, const char* alignment_filepath,
                               int64_t (*kmer_to_dp_id_func) (char*), bool has_header) {
    
    stList* signal_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);
    
    FILE* align_file = fopen(alignment_filepath, "r");
    
    stList* tokens;
    int64_t line_length;
    char* kmer;
    char* signal_str;
    int64_t* dp_id_ptr;
    double* signal_ptr;
    bool warned = false;
    
    char* line = stFile_getLineFromFile(align_file);
    if (has_header) {
        line = stFile_getLineFromFile(align_file);
    }
    while (line != NULL) {
        tokens = stString_split(line);
        line_length = stList_length(tokens);
        
        if (!warned) {
            if (line_length != NUM_ALIGNMENT_COLS) {
                fprintf(stderr, "Input format has changed from design period, HDP may receive incorrect data.");
                warned = true;
            }
        }
        
        signal_str = (char*) stList_get(tokens, ALIGNMENT_SIGNAL_COL);
        kmer = (char*) stList_get(tokens, ALIGNMENT_KMER_COL);
        
        signal_ptr = (double*) malloc(sizeof(double));
        dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));
        
        sscanf(signal_str, "%lf", signal_ptr);
        *dp_id_ptr = kmer_to_dp_id_func(kmer);
        
        stList_append(signal_list, signal_ptr);
        stList_append(dp_id_list, dp_id_ptr);
        
        stList_destruct(tokens);
        free(line);
        line = stFile_getLineFromFile(align_file);
    }
    
    fclose(align_file);
    
    int64_t data_length;
    
    double* signal = stList_toDoublePtr(signal_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &data_length);
    
    stList_destruct(signal_list);
    stList_destruct(dp_id_list);
    
    reset_hdp_data(hdp);
    pass_data_to_hdp(hdp, signal, dp_ids, data_length);
}



// n^k
int64_t power(int64_t n, int64_t k) {
    int64_t num = 1;
    
    for (int64_t i = 0; i < k; i++) {
        num *= n;
    }
    
    return num;
}

//  ((n k))
int64_t multiset_number(int64_t n, int64_t k) {
    int64_t num = 1;
    for (int64_t m = n + k - 1; m >= n; m--) {
        num *= m;
    }
    for (int64_t m = k; m >= 2; m--) {
        num /= m;
    }
    return num;
}

int64_t* get_word(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* word = (int64_t*) malloc(sizeof(int64_t) * word_length);
    int64_t id_remainder = word_id;
    for (int64_t i = 0; i < word_length; i++) {
        word[word_length - i - 1] = id_remainder % alphabet_size;
        id_remainder /= alphabet_size;
    }
    return word;
}

int64_t* get_word_multiset(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* multiset = get_word(word_id, alphabet_size, word_length);
    
    // selection sort 'cause whatever
    int64_t min_idx;
    int64_t temp;
    for (int64_t i = 0; i < word_length; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < word_length; j++) {
            if (multiset[j] < multiset[min_idx]) {
                min_idx = j;
            }
        }
        temp = multiset[i];
        multiset[i] = multiset[min_idx];
        multiset[min_idx] = temp;
    }
    
    return multiset;
}

int64_t multiset_id_internal(int64_t* tail, int64_t tail_length, int64_t alphabet_min, int64_t alphabet_size) {
    int64_t head = tail[0];
    if (tail_length == 1) {
        return head - alphabet_min;
    }
    int64_t step = 0;
    for (int64_t i = alphabet_min; i < alphabet_size; i++) {
        if (head > i) {
            step += multiset_number(alphabet_size - i, tail_length - 1);
        }
        else {
            return step + multiset_id_internal(&(tail[1]), tail_length - 1, i, alphabet_size);
        }
    }
    fprintf(stderr, "Character outside alphabet included in multiset\n");
    exit(EXIT_FAILURE);
}

int64_t multiset_id(int64_t* multiset, int64_t length, int64_t alphabet_size) {
    return multiset_id_internal(multiset, length, 0, alphabet_size);
}

int64_t word_id_to_multiset_id(int64_t word_id, int64_t alphabet_size, int64_t word_length) {
    int64_t* multiset = get_word_multiset(word_id, alphabet_size, word_length);
    int64_t id = multiset_id(multiset, word_length, alphabet_size);
    free(multiset);
    return id;
}


                                        
                                        
                                        
                                        

int64_t flat_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    return num_leaves + 1;
}

void flat_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    int64_t last_dp_id = power(alphabet_size, kmer_length);
    
    for (int64_t id = 0; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}

HierarchicalDirichletProcess* flat_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                             double leaf_gamma, double sampling_grid_start,
                                             double sampling_grid_stop, int64_t sampling_grid_length,
                                             const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 2);
    gamma_params[0] = base_gamma;
    gamma_params[1] = leaf_gamma;
    
    int64_t num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 2, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* flat_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                               double base_gamma_alpha, double base_gamma_beta,
                                               double leaf_gamma_alpha, double leaf_gamma_beta,
                                               double sampling_grid_start, double sampling_grid_stop,
                                               int64_t sampling_grid_length,
                                               const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 2);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = leaf_gamma_alpha;
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 2);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = leaf_gamma_beta;
    
    int64_t num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 2, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

int64_t multiset_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(kmer_length, alphabet_size);
    return num_leaves + num_middle_dps + 1;
}

void multiset_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = multiset_number(kmer_length, alphabet_size);
    
    // set kmer parents to multisets
    int64_t multiset_id;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        multiset_id = word_id_to_multiset_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + multiset_id);
    }
    
    // set multiset parents to base dp
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

HierarchicalDirichletProcess* multiset_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                                 double middle_gamma, double leaf_gamma,
                                                 double sampling_grid_start,
                                                 double sampling_grid_stop, int64_t sampling_grid_length,
                                                 const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* multiset_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                                   double base_gamma_alpha, double base_gamma_beta,
                                                   double middle_gamma_alpha, double middle_gamma_beta,
                                                   double leaf_gamma_alpha, double leaf_gamma_beta,
                                                   double sampling_grid_start, double sampling_grid_stop,
                                                   int64_t sampling_grid_length,
                                                   const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

int64_t middle_2_nts_hdp_num_dps(int64_t alphabet_size, int64_t kmer_length) {
    if (kmer_length <= 2) {
        fprintf(stderr, "k-mer is not long enough for middle 2 nucleotides HDP\n");
        exit(EXIT_FAILURE);
    }
    return power(alphabet_size, kmer_length) + power(alphabet_size, 2) + 1;
}

int64_t kmer_id_to_middle_nts_id(int64_t kmer_id, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* kmer = get_word(kmer_id, alphabet_size, kmer_length);
    int64_t id = alphabet_size * kmer[kmer_length / 2 - 1] + kmer[kmer_length / 2];
    free(kmer);
    return id;
}

void middle_2_nts_hdp_model_internal(HierarchicalDirichletProcess* hdp, int64_t alphabet_size, int64_t kmer_length) {
    
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = power(alphabet_size, 2);
    
    int64_t middle_dp_id;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        middle_dp_id = kmer_id_to_middle_nts_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, middle_dp_id);
    }
    
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t id = num_leaves; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}

HierarchicalDirichletProcess* middle_2_nts_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                                     double middle_gamma, double leaf_gamma,
                                                     double sampling_grid_start,
                                                     double sampling_grid_stop, int64_t sampling_grid_length,
                                                     const char* model_filepath) {
    if (kmer_length % 2 != 0) {
        fprintf(stderr, "Warning: middle 2 nucleotides of odd length kmer is ambiguous. Resolving arbitrarily.\n");
    }
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* middle_2_nts_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                                       double base_gamma_alpha, double base_gamma_beta,
                                                       double middle_gamma_alpha, double middle_gamma_beta,
                                                       double leaf_gamma_alpha, double leaf_gamma_beta,
                                                       double sampling_grid_start, double sampling_grid_stop,
                                                       int64_t sampling_grid_length,
                                                       const char* model_filepath) {
    if (kmer_length % 2 != 0) {
        fprintf(stderr, "Warning: middle 2 nucleotides of odd length kmer is ambiguous. Resolving arbitrarily.\n");
    }
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}
