//
//  nanopore_hdp.c
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

#include <stdio.h>
#include <stdbool.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "nanopore_hdp.h"
#include "sonLib.h"

// in 0-based index
#declare KMER_COL 9
#declare SIGNAL_COL 13
#declare NUM_COLS 14

void normal_inverse_gamma_params_from_minION(const char* signal_lookup_table, double* mu_out, double* nu_out,
                                             double* alpha_out, double* beta_out) {
    
    FILE* lookup_table = fopen(signal_lookup_table, "r");
    
    double* means;
    int length;
    
    //TODO: get the mean vector
    
    normal_inverse_gamma_params(means, length, mu_out, nu_out, alpha_out, beta_out);
    
}

// fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int num_dps, int depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int sampling_grid_length,
                                         const char* signal_lookup_table_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(signal_lookup_table_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc(num_dps, depth, gamma, sampling_grid_start, sampling_grid_stop,
                             sampling_grid_length, mu, nu, alpha, beta);
}

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int num_dps, int depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int sampling_grid_length,
                                           const char* signal_lookup_table_filepath) {
    
    double mu, nu, alpha, beta;
    normal_inverse_gamma_params_from_minION(signal_lookup_table_filepath, &mu, &nu, &alpha, &beta);
    return new_hier_dir_proc_2(num_dps, depth, gamma_alpha, gamma_beta, sampling_grid_start,
                               sampling_grid_stop, sampling_grid_length, mu, nu, alpha, beta);
}

int* stList_toIntPtr(stList* list, int* length_out) {
    int length = (int) stList_length(list);
    int* int_arr = (int*) malloc(sizeof(int) * length);
    int* entry;
    for (int i = 0; i < length; i++) {
        entry = (int*) stList_get(list, i);
        int_arr[0] = *entry;
    }
    *length_out = length;
    return int_arr;
}

double* stList_toDoublePtr(stList* list, int* length_out) {
    int length  = stList_length(list);
    double* double_arr = (double*) malloc(sizeof(double) * length);
    double* entry;
    for (int i = 0; i < length; i++) {
        entry = (double*) stList_get(list, i);
        double_arr[0] = *entry;
    }
    *length_out = length;
    return double_arr;
}

void update_hdp_from_alignment(HierarchicalDirichletProcess* hdp, const char* alignment_filepath,
                               int (*kmer_to_dp_id_func) (char*), bool has_header) {
    
    stList* signal_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);
    
    FILE* align_file = fopen(alignment_filepath, "r");
    
    stList* tokens;
    int64_t line_length;
    char* kmer;
    char* signal_str;
    int* dp_id_ptr;
    double* signal_ptr;
    bool warned = false;
    
    char* line = stFile_getLineFromFile(align_file);
    if (has_header) {
        line = stFile_getLineFromFile(align_file);
    }
    while (line != NULL) {
        tokens = stString_split(line);
        line_length = stList_getLength(tokens);
        
        if (!warned) {
            if (line_length != NUM_COLS) {
                fprintf(stderr, "Input format has changed from design period, HDP may receive incorrect data.");
                warned = true;
            }
        }
        
        signal_str = (char*) stList_get(tokens, SIGNAL_COL);
        kmer = (char*) stList_get(tokens, KMER_COL);
        
        signal_ptr = (double*) malloc(sizeof(double));
        dp_id_ptr = (double*) malloc(sizeof(int));
        
        sscanf(signal_str, "%lf", signal_ptr);
        *dp_id_ptr = kmer_to_dp_id_func(kmer);
        
        stList_append(signal_list, signal_ptr);
        stList_append(dp_id_list, dp_id_ptr);
        
        stList_destruct(tokens);
        free(line);
        line = stFile_getLineFromFile(align_file);
    }
    
    fclose(align_file);
    
    int data_length;
    
    double* signal = stList_toDoublePtr(signal_list, &data_length);
    int* dp_ids = stList_toIntPtr(dp_id_list, &data_length);
    
    stList_destruct(signal_list);
    stList_destruct(dp_id_list);
    
    reset_hdp_data(hdp);
    pass_data_to_hdp(hdp, signal, dp_ids, data_length);
}

// n^k
int power(int n, int k) {
    int num = 1;
    
    for (int i = 0; i < k; i++) {
        num *= n;
    }
    
    return num;
}

//  ((n k))
int multiset_number(int n, int k) {
    int num = 1;
    for (int m = n + k - 1; m >= n; m--) {
        num *= m;
    }
    for (int m = k; m >= 2; m--) {
        num /= m;
    }
    return num;
}

int flat_hdp_num_dps(int alphabet_size, int kmer_length) {
    int num_leaves = power(alphabet_size, kmer_length);
    return = num_leaves + 1;
}

void flat_hdp_model_internal(HierarchicalDirichletProcess* hdp) {
    int last_dp_id = hdp->num_dps - 1;
    
    for (int id = 0; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}


HierarchicalDirichletProcess* flat_hdp_model(int alphabet_size, int kmer_length, double base_gamma,
                                             double leaf_gamma, double sampling_grid_start,
                                             double sampling_grid_stop, int sampling_grid_length,
                                             const char* signal_lookup_table_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 2);
    gamma_params[0] = base_gamma;
    gamma_params[1] = leaf_gamma;
    
    int num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 2, gamma, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   signal_lookup_table_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* flat_hdp_model_2(int alphabet_size, int kmer_length,
                                               double base_gamma_alpha, double base_gamma_beta,
                                               double leaf_gamma_alpha, double leaf_gamma_beta,
                                               double sampling_grid_start, double sampling_grid_stop,
                                               int sampling_grid_length,
                                               const char* signal_lookup_table_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 2);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = leaf_gamma_alpha;
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 2);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = leaf_gamma_beta;
    
    int num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 2, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     signal_lookup_table_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

int multiset_hdp_num_dps(int alphabet_size, int kmer_length) {
    int num_leaves = power(alphabet_size, kmer_length);
    int num_middle_dps = multiset_number(kmer_length, alphabet_size);
    return num_leaves + num_middle_dps + 1;
}

int* get_kmer(int kmer_id, int alphabet_size, int kmer_length) {
    int* kmer = (int*) malloc(sizeof(int) * kmer_length);
    int remainder = kmer_id;
    for (int i = 0; i < kmer_length; i++) {
        kmer[kmer_length - i - 1] = remainder % alphabet_size;
        remainder /= alphabet_size;
    }
    return kmer;
}

int* get_kmer_multiset(int kmer_id, int alphabet_size, int kmer_length) {
    int* multiset = get_kmer(kmer_id, alphabet_size, kmer_length);
    
    // selection sort 'cause whatever
    int min_idx;
    int temp;
    for (int i = 0; i < kmer_length; i++) {
        min_idx = i;
        for (int j = i + 1; j < kmer_length; j++) {
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

int multiset_id(int* tail, int tail_length, int alphabet_min, int alphabet_size) {
    int head = tail[0];
    if (tail_length == 1) {
        return head - alphabet_min;
    }
    int step = 0;
    for (int i = alphabet_min; i < alphabet_size; i++) {
        if (head > i) {
            step += multiset_number(alphabet_size - i, tail_length - 1);
        }
        else {
            return step + multiset_id_internal(&tail[1], tail_length - 1, i, alphabet_size);
        }
    }
    fprintf(stderr, "Character outside alphabet included in multiset\n");
    exit(EXIT_FAILURE);
}



int kmer_id_to_multiset_id(int kmer_id, int alphabet_size, int kmer_length) {
    int* multiset = get_kmer_multiset(kmer_id, alphabet_size, kmer_length);
    int id = multiset_id(multiset, kmer_length, 0, alphabet_size);
    free(multiset);
    return id;
}

void multiset_hdp_model_internal(HierarchicalDirichletProcess* hdp, int alphabet_size, int kmer_length) {
    int num_leaves = power(alphabet_size, kmer_length);
    
    // set kmer parents to multisets
    int multiset_id;
    for (int kmer_id = 0; kmer_id < num_leaves; id++) {
        multiset_id = kmer_id_to_multiset_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + multiset_id);
    }
    
    // set multiset parents to base dp
    int last_dp_id = hdp->num_dps - 1;
    for (int middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

HierarchicalDirichletProcess* multiset_hdp_model(int alphabet_size, int kmer_length, double base_gamma,
                                                 double middle_gamma, double leaf_gamma,
                                                 double sampling_grid_start,
                                                 double sampling_grid_stop, int sampling_grid_length,
                                                 const char* signal_lookup_table_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   signal_lookup_table_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* multiset_hdp_model_2(int alphabet_size, int kmer_length,
                                                   double base_gamma_alpha, double base_gamma_beta,
                                                   double middle_gamma_alpha, double middle_gamma_beta,
                                                   double leaf_gamma_alpha, double leaf_gamma_beta,
                                                   double sampling_grid_start, double sampling_grid_stop,
                                                   int sampling_grid_length,
                                                   const char* signal_lookup_table_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int num_dps = multiset_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     signal_lookup_table_filepath);
    
    multiset_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

int middle_2_nts_hdp_num_dps(alphabet_size, kmer_length) {
    if (kmer_length <= 2) {
        fprintf(stderr, "k-mer is not long enough for middle 2 nucleotides HDP\n");
        exit(EXIT_FAILURE);
    }
    
    return power(alphabet_size, kmer_length) + power(alphabet_size, 2) + 1;
}

int kmer_id_to_middle_nts_id(int kmer_id, int alphabet_size, int kmer_length) {
    int* kmer = get_kmer(kmer_id, alphabet_size, kmer_length);
    int id = alphabet_size * kmer[kmer_length / 2 - 1] + kmer[kmer_length / 2];
    free(kmer);
    return id;
}

void middle_2_nts_hdp_model_internal(HierarchicalDirichletProcess* hdp, int alphabet_size, int kmer_length) {
    
    int num_leaves = power(alphabet_size, kmer_length);
    
    int middle_dp_id;
    for (int kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        middle_dp_id = kmer_id_to_middle_nts_id(kmer_id, alphabet_size, kmer_length);
        set_dir_proc_parent(hdp, kmer_id, middle_dp_id)
    }
    
    int last_dp_id = hdp->num_dps - 1;
    for (int id = num_leaves; id < last_dp_id; id++) {
        set_dir_proc_parent(hdp, id, last_dp_id);
    }
}

HierarchicalDirichletProcess* middle_2_nts_hdp_model(int alphabet_size, int kmer_length, double base_gamma,
                                                     double middle_gamma, double leaf_gamma,
                                                     double sampling_grid_start,
                                                     double sampling_grid_stop, int sampling_grid_length,
                                                     const char* signal_lookup_table_filepath) {
    if (kmer_length %2 != 2) {
        fprintf(stderr, "Warning: middle 2 nucleotides of odd length kmer is ambiguous. Resolving arbitrarily.\n");
    }
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   signal_lookup_table_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}

HierarchicalDirichletProcess* middle_2_nts_hdp_model(int alphabet_size, int kmer_length,
                                                     double base_gamma_alpha, double base_gamma_beta,
                                                     double middle_gamma_alpha, double middle_gamma_beta,
                                                     double leaf_gamma_alpha, double leaf_gamma_beta,
                                                     double sampling_grid_start, double sampling_grid_stop,
                                                     int sampling_grid_length,
                                                     const char* signal_lookup_table_filepath) {
    if (kmer_length %2 != 2) {
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
    
    int num_dps = middle_2_nts_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     signal_lookup_table_filepath);
    
    middle_2_nts_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    return hdp;
}



