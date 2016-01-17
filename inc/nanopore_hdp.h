//
//  nanopore_hdp.h
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

#ifndef nanopore_hdp_h
#define nanopore_hdp_h

#include <stdbool.h>
#include <inttypes.h>
#include "hdp.h"

// for fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int64_t sampling_grid_length,
                                         const char* model_filepath);

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int64_t num_dps, int64_t depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int64_t sampling_grid_length,
                                           const char* model_filepath);

// single level HDP
HierarchicalDirichletProcess* flat_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                             double leaf_gamma, double sampling_grid_start,
                                             double sampling_grid_stop, int64_t sampling_grid_length,
                                             const char* model_filepath);
HierarchicalDirichletProcess* flat_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                              double base_gamma_alpha, double base_gamma_beta,
                                              double leaf_gamma_alpha, double leaf_gamma_beta,
                                              double sampling_grid_start, double sampling_grid_stop,
                                              int64_t sampling_grid_length,
                                              const char* model_filepath);

// second level of HDP based on multiset of nucleotides
HierarchicalDirichletProcess* multiset_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                                 double middle_gamma, double leaf_gamma,
                                                 double sampling_grid_start,
                                                 double sampling_grid_stop, int64_t sampling_grid_length,
                                                 const char* model_filepath);
HierarchicalDirichletProcess* multiset_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                                   double base_gamma_alpha, double base_gamma_beta,
                                                   double middle_gamma_alpha, double middle_gamma_beta,
                                                   double leaf_gamma_alpha, double leaf_gamma_beta,
                                                   double sampling_grid_start, double sampling_grid_stop,
                                                   int64_t sampling_grid_length,
                                                   const char* model_filepath);

// second level of HDP based on middle 2 nucleotides
HierarchicalDirichletProcess* middle_2_nts_hdp_model(int64_t alphabet_size, int64_t kmer_length, double base_gamma,
                                                     double middle_gamma, double leaf_gamma,
                                                     double sampling_grid_start,
                                                     double sampling_grid_stop, int64_t sampling_grid_length,
                                                     const char* model_filepath);
HierarchicalDirichletProcess* middle_2_nts_hdp_model_2(int64_t alphabet_size, int64_t kmer_length,
                                                       double base_gamma_alpha, double base_gamma_beta,
                                                       double middle_gamma_alpha, double middle_gamma_beta,
                                                       double leaf_gamma_alpha, double leaf_gamma_beta,
                                                       double sampling_grid_start, double sampling_grid_stop,
                                                       int64_t sampling_grid_length,
                                                       const char* model_filepath);

// you write your own function that maps kmers to integers
void update_hdp_from_alignment(HierarchicalDirichletProcess* hdp, const char* alignment_filepath,
                               int64_t (*kmer_to_dp_id_func) (char*), bool has_header);


// n^k
int64_t power(int64_t n, int64_t k);
//  ((n k))
int64_t multiset_number(int64_t n, int64_t k);
// get word as int array by lexicographic order index
int64_t* get_word(int64_t word_id, int64_t alphabet_size, int64_t word_length);
// get multiset of a word as sorted int array by lexicographic order indx
int64_t* get_word_multiset(int64_t word_id, int64_t alphabet_size, int64_t word_length);
// get multiset lexicographic index from a multiset as sorted int array
int64_t multiset_id(int64_t* multiset, int64_t length, int64_t alphabet_size);
// get lexicographic index of multiset from lexicographic index of word
int64_t word_id_to_multiset_id(int64_t word_id, int64_t alphabet_size, int64_t word_length);

#endif /* nanopore_hdp_h */
