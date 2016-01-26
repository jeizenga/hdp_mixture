//
//  nanopore_hdp.c
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

// in 0-based index
#define ALIGNMENT_KMER_COL 9
#define ALIGNMENT_STRAND_COL 4
#define ALIGNMENT_SIGNAL_COL 13
#define NUM_ALIGNMENT_COLS 15

#define MODEL_ROW_HEADER_LENGTH 1
#define MODEL_MEAN_ENTRY 0
#define MODEL_ENTRY_LENGTH 5

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "hdp.h"
#include "hdp_math_utils.h"
#include "nanopore_hdp.h"
#include "sonLib.h"

struct NanoporeHDP {
    HierarchicalDirichletProcess* hdp;
    char* alphabet;
    int64_t alphabet_size;
    int64_t kmer_length;
};

NanoporeHDP* package_nanopore_hdp(HierarchicalDirichletProcess* hdp, const char* alphabet, int64_t alphabet_size,
                                  int64_t kmer_length) {
    
    NanoporeHDP* nhdp = (NanoporeHDP*) malloc(sizeof(NanoporeHDP));
    
    // copy and sort alphabet
    char* internal_alphabet = (char*) malloc(sizeof(char) * (alphabet_size + 1));
    for (int64_t i = 0; i < alphabet_size; i++) {
        internal_alphabet[i] = alphabet[i];
    }
    
    int64_t min_idx;
    char temp;
    for (int64_t i = 0; i < alphabet_size; i++) {
        min_idx = i;
        for (int64_t j = i + 1; j < alphabet_size; j++) {
            if (internal_alphabet[j] < internal_alphabet[min_idx]) {
                min_idx = j;
            }
        }
        temp = internal_alphabet[i];
        internal_alphabet[i] = internal_alphabet[min_idx];
        internal_alphabet[min_idx] = temp;
    }
    
    for (int64_t i = 1; i < alphabet_size; i++) {
        if (alphabet[i - 1] == alphabet[i]) {
            fprintf(stderr, "Characters of alphabet must be distinct.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    internal_alphabet[alphabet_size] = '\0';
    
    nhdp->hdp = hdp;
    nhdp->alphabet = internal_alphabet;
    nhdp->alphabet_size = alphabet_size;
    nhdp->kmer_length = kmer_length;
    
    return nhdp;
}

void destroy_nanopore_hdp(NanoporeHDP* nhdp) {
    destroy_hier_dir_proc(nhdp->hdp);
    free(nhdp->alphabet);
    free(nhdp);
}

int64_t get_nanopore_hdp_kmer_length(NanoporeHDP* nhdp) {
    return nhdp->kmer_length;
}


int64_t get_nanopore_hdp_alphabet_size(NanoporeHDP* nhdp) {
    return nhdp->alphabet_size;
}

char* get_nanopore_hdp_alphabet(NanoporeHDP* nhdp) {
    char* alphabet = nhdp->alphabet;
    int64_t alphabet_size = nhdp->alphabet_size;
    char* copy = (char*) malloc(sizeof(char) * (alphabet_size + 1));
    for (int64_t i = 0; i < alphabet_size; i++) {
        copy[i] = alphabet[i];
    }
    copy[alphabet_size] = '\0';
    return copy;
}


// wrappers
void execute_nhdp_gibbs_sampling(NanoporeHDP* nhdp, int64_t num_samples, int64_t burn_in,
                                 int64_t thinning, bool verbose) {
    execute_gibbs_sampling(nhdp->hdp, num_samples, burn_in, thinning, verbose);
}

void execute_nhdp_gibbs_sampling_with_snapshots(NanoporeHDP* nhdp,
                                                int64_t num_samples, int64_t burn_in, int64_t thinning,
                                                void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                                void* snapshot_func_args, bool verbose) {
    execute_gibbs_sampling_with_snapshots(nhdp->hdp, num_samples, burn_in, thinning, snapshot_func, snapshot_func_args,
                                          verbose);
}
void finalize_nhdp_distributions(NanoporeHDP* nhdp) {
    finalize_distributions(nhdp->hdp);
}

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

void update_nhdp_from_alignment(NanoporeHDP* nhdp, const char* alignment_filepath, bool has_header) {
    update_nhdp_from_alignment_with_filter(nhdp, alignment_filepath, has_header, NULL);
}

void update_nhdp_from_alignment_with_filter(NanoporeHDP* nhdp, const char* alignment_filepath,
                                            bool has_header, const char* strand_filter) {
    
    stList* signal_list = stList_construct3(0, &free);
    stList* dp_id_list = stList_construct3(0, &free);
    
    FILE* align_file = fopen(alignment_filepath, "r");
    if (align_file == NULL) {
        fprintf(stderr, "Alignment %s file does not exist.\n", alignment_filepath);
        exit(EXIT_FAILURE);
    }
    
    stList* tokens;
    int64_t line_length;
    char* kmer;
    char* strand;
    char* signal_str;
    int64_t* dp_id_ptr;
    double* signal_ptr;
    bool warned = false;
    int proceed = 0;
    
    char* line = stFile_getLineFromFile(align_file);
    if (has_header) {
        line = stFile_getLineFromFile(align_file);
    }
    while (line != NULL) {
        tokens = stString_split(line);
        line_length = stList_length(tokens);
        
        if (!warned) {
            if (line_length != NUM_ALIGNMENT_COLS) {
                fprintf(stderr, "Input format has changed from design period, HDP may receive incorrect data.\n");
                warned = true;
            }
        }
        
        strand = (char*) stList_get(tokens, ALIGNMENT_STRAND_COL);
        
        if (strand_filter != NULL) {
            proceed = strcmp(strand, strand_filter);
        }
        
        if (proceed == 0) {
            signal_str = (char*) stList_get(tokens, ALIGNMENT_SIGNAL_COL);
            kmer = (char*) stList_get(tokens, ALIGNMENT_KMER_COL);
            
            signal_ptr = (double*) malloc(sizeof(double));
            dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));
            sscanf(signal_str, "%lf", signal_ptr);
            *dp_id_ptr = kmer_id(kmer, nhdp->alphabet, nhdp->alphabet_size, nhdp->kmer_length);
            
            stList_append(signal_list, signal_ptr);
            stList_append(dp_id_list, dp_id_ptr);
        }
        
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
    
    reset_hdp_data(nhdp->hdp);
    pass_data_to_hdp(nhdp->hdp, signal, dp_ids, data_length);
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

int64_t word_id(int64_t* word, int64_t alphabet_size, int64_t word_length) {
    int64_t id = 0;
    int64_t step = 1;
    for (int64_t i = word_length - 1; i >= 0; i--) {
        id += step * word[i];
        step *= alphabet_size;
    }
    return id;
}

int64_t* kmer_to_word(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* word = (int64_t*) malloc(sizeof(int64_t) * kmer_length);
    for (int64_t i = 0; i < kmer_length; i++) {
        int64_t j = 0;
        while (kmer[i] != alphabet[j]) {
            j++;
            if (j == alphabet_size) {
                fprintf(stderr, "K-mer contains character outside alphabet.\n");
                exit(EXIT_FAILURE);
            }
        }
        word[i] = j;
    }
    return word;
}

int64_t kmer_id(char* kmer, char* alphabet, int64_t alphabet_size, int64_t kmer_length) {
    int64_t* word = kmer_to_word(kmer, alphabet, alphabet_size, kmer_length);
    int64_t id = word_id(word, alphabet_size, kmer_length);
    free(word);
    return id;
}

int64_t standard_kmer_id(char* kmer, int64_t kmer_length) {
    return kmer_id(kmer, "ACGT", 4, kmer_length);
}

double get_nanopore_kmer_density(NanoporeHDP* nhdp, double x, char* kmer) {
    int64_t dp_id = kmer_id(kmer, nhdp->alphabet, nhdp->alphabet_size, nhdp->kmer_length);
    return dir_proc_density(nhdp->hdp, x, dp_id);
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

NanoporeHDP* flat_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                            double leaf_gamma, double sampling_grid_start, double sampling_grid_stop,
                            int64_t sampling_grid_length, const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 2);
    gamma_params[0] = base_gamma;
    gamma_params[1] = leaf_gamma;
    
    int64_t num_dps = flat_hdp_num_dps(alphabet_size, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 2, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    flat_hdp_model_internal(hdp, alphabet_size, kmer_length);
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

NanoporeHDP* flat_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                              double base_gamma_alpha, double base_gamma_beta, double leaf_gamma_alpha,
                              double leaf_gamma_beta, double sampling_grid_start, double sampling_grid_stop,
                              int64_t sampling_grid_length, const char* model_filepath) {
    
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
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
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

NanoporeHDP* multiset_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                double base_gamma, double middle_gamma, double leaf_gamma,
                                double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
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
    
    finalize_hdp_structure(hdp);

    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

NanoporeHDP* multiset_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                  double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                  double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                  double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
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
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
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

NanoporeHDP* middle_2_nts_hdp_model(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                    double base_gamma, double middle_gamma, double leaf_gamma,
                                    double sampling_grid_start, double sampling_grid_stop, int64_t sampling_grid_length,
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
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

NanoporeHDP* middle_2_nts_hdp_model_2(const char* alphabet, int64_t alphabet_size, int64_t kmer_length,
                                      double base_gamma_alpha, double base_gamma_beta, double middle_gamma_alpha,
                                      double middle_gamma_beta, double leaf_gamma_alpha, double leaf_gamma_beta,
                                      double sampling_grid_start, double sampling_grid_stop,
                                      int64_t sampling_grid_length, const char* model_filepath) {
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
    
    finalize_hdp_structure(hdp);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    return nhdp;
}

int64_t purine_composition_hdp_num_dps(int64_t num_purines, int64_t num_pyrimidines, int64_t kmer_length) {
    int64_t num_leaves = power(num_purines + num_pyrimidines, kmer_length);
    int64_t num_middle_dps = kmer_length + 1;
    return num_leaves + num_middle_dps + 1;
}

void purine_composition_hdp_model_internal(HierarchicalDirichletProcess* hdp, bool* purine_alphabet,
                                           int64_t alphabet_size, int64_t kmer_length) {
    int64_t num_leaves = power(alphabet_size, kmer_length);
    int64_t num_middle_dps = kmer_length + 1;
    
    // set kmer parents to purine multisets
    int64_t num_purines;
    int64_t* word;
    for (int64_t kmer_id = 0; kmer_id < num_leaves; kmer_id++) {
        word = get_word(kmer_id, alphabet_size, kmer_length);
        num_purines = 0;
        for (int64_t i = 0; i < kmer_length; i++) {
            if (purine_alphabet[word[i]]) {
                num_purines++;
            }
        }
        free(word);
        set_dir_proc_parent(hdp, kmer_id, num_leaves + num_purines);
    }
    
    // set purine set parents to base dp
    int64_t last_dp_id = num_leaves + num_middle_dps;
    for (int64_t middle_dp_id = num_leaves; middle_dp_id < last_dp_id; middle_dp_id++) {
        set_dir_proc_parent(hdp, middle_dp_id, last_dp_id);
    }
}

NanoporeHDP* purine_composition_hdp_model(char* purine_alphabet, int64_t num_purines,
                                          char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                          int64_t kmer_length, double base_gamma, double middle_gamma,
                                          double leaf_gamma, double sampling_grid_start, double sampling_grid_stop,
                                          int64_t sampling_grid_length, const char* model_filepath) {
    
    double* gamma_params = (double*) malloc(sizeof(double) * 3);
    gamma_params[0] = base_gamma;
    gamma_params[1] = middle_gamma;
    gamma_params[2] = leaf_gamma;
    
    int64_t num_dps = purine_composition_hdp_num_dps(num_purines, num_pyrimidines, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp(num_dps, 3, gamma_params, sampling_grid_start,
                                                   sampling_grid_stop, sampling_grid_length,
                                                   model_filepath);
    
    int64_t alphabet_size = num_purines + num_pyrimidines;
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        alphabet[i] = purine_alphabet[i];
    }
    for (int64_t i = 0; i < num_pyrimidines; i++) {
        alphabet[i + num_purines] = pyrimidine_alphabet[i];
    }
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    // get back the alphabet in the internal ordering
    free(alphabet);
    alphabet = get_nanopore_hdp_alphabet(nhdp);
    bool* purines = (bool*) malloc(sizeof(bool) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        purine_alphabet[i] = false;
        for (int64_t j = 0; j < num_purines; j++) {
            if (alphabet[i] == purine_alphabet[j]) {
                purines[i] = true;
                break;
            }
        }
    }
    free(alphabet);
    
    purine_composition_hdp_model_internal(hdp, purines, alphabet_size, kmer_length);
    free(purines);
    
    finalize_hdp_structure(hdp);
    
    return nhdp;
}

NanoporeHDP* purine_composition_hdp_model_2(char* purine_alphabet, int64_t num_purines,
                                            char* pyrimidine_alphabet, int64_t num_pyrimidines,
                                            int64_t kmer_length, double base_gamma_alpha, double base_gamma_beta,
                                            double middle_gamma_alpha, double middle_gamma_beta,
                                            double leaf_gamma_alpha, double leaf_gamma_beta, double sampling_grid_start,
                                            double sampling_grid_stop, int64_t sampling_grid_length,
                                            const char* model_filepath) {
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * 3);
    gamma_alpha[0] = base_gamma_alpha;
    gamma_alpha[1] = middle_gamma_alpha;
    gamma_alpha[2] = leaf_gamma_alpha;
    
    
    double* gamma_beta = (double*) malloc(sizeof(double) * 3);
    gamma_beta[0] = base_gamma_beta;
    gamma_beta[1] = middle_gamma_beta;
    gamma_beta[2] = leaf_gamma_beta;
    
    int64_t num_dps = purine_composition_hdp_num_dps(num_purines, num_pyrimidines, kmer_length);
    
    HierarchicalDirichletProcess* hdp = minION_hdp_2(num_dps, 3, gamma_alpha, gamma_beta, sampling_grid_start,
                                                     sampling_grid_stop, sampling_grid_length,
                                                     model_filepath);
    
    int64_t alphabet_size = num_purines + num_pyrimidines;
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    for (int64_t i = 0; i < num_purines; i++) {
        alphabet[i] = purine_alphabet[i];
    }
    for (int64_t i = 0; i < num_pyrimidines; i++) {
        alphabet[i + num_purines] = pyrimidine_alphabet[i];
    }
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
    
    // get back the alphabet in the internal ordering
    free(alphabet);
    alphabet = get_nanopore_hdp_alphabet(nhdp);
    bool* purines = (bool*) malloc(sizeof(bool) * alphabet_size);
    for (int64_t i = 0; i < alphabet_size; i++) {
        purines[i] = false;
        for (int64_t j = 0; j < num_purines; j++) {
            if (alphabet[i] == purine_alphabet[j]) {
                purines[i] = true;
                break;
            }
        }
    }
    
    free(alphabet);
    
    purine_composition_hdp_model_internal(hdp, purines, alphabet_size, kmer_length);
    free(purines);
    
    finalize_hdp_structure(hdp);
    
    return nhdp;
}

void serialize_nhdp(NanoporeHDP* nhdp, const char* filepath) {
    FILE* out = fopen(filepath, "w");
    
    fprintf(out, "%"PRId64"\n", nhdp->alphabet_size);
    fprintf(out, "%s\n", nhdp->alphabet);
    fprintf(out, "%"PRId64"\n", nhdp->kmer_length);
    serialize_hdp(nhdp->hdp, out);
    
    fclose(out);
}

NanoporeHDP* deserialize_nhdp(const char* filepath) {
    FILE* in = fopen(filepath, "r");
    
    char* line = stFile_getLineFromFile(in);
    int64_t alphabet_size;
    sscanf(line, "%"SCNd64, &alphabet_size);
    free(line);
    
    line = stFile_getLineFromFile(in);
    char* alphabet = (char*) malloc(sizeof(char) * alphabet_size);
    sscanf(line, "%s", alphabet);
    free(line);
    
    line = stFile_getLineFromFile(in);
    int64_t kmer_length;
    sscanf(line, "%"SCNd64, &kmer_length);
    free(line);
    
    HierarchicalDirichletProcess* hdp = deserialize_hdp(in);
    
    fclose(in);
    
    NanoporeHDP* nhdp = package_nanopore_hdp(hdp, alphabet, alphabet_size, kmer_length);
                                    
    free(alphabet);
    
    return nhdp;
}


