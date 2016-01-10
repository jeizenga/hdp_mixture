//
//  nanopore_hdp.h
//  
//
//  Created by Jordan Eizenga on 1/8/16.
//
//

#ifndef nanopore_hdp_h
#define nanopore_hdp_h

// for fixed concentration parameters 'gamma' for each depth
HierarchicalDirichletProcess* minION_hdp(int num_dps, int depth, double* gamma, double sampling_grid_start,
                                         double sampling_grid_stop, int sampling_grid_length,
                                         const char* signal_lookup_table_filepath);

// Gamma distribution prior on the concentration parameters 'gamma'
// must designate vector of 'alpha' and 'beta' parameters of distribution for each depth
HierarchicalDirichletProcess* minION_hdp_2(int num_dps, int depth, double* gamma_alpha,
                                           double* gamma_beta, double sampling_grid_start,
                                           double sampling_grid_stop, int sampling_grid_length,
                                           const char* signal_lookup_table_filepath);

// you write your own function that maps kmers to integers
void update_hdp_from_alignment(HierarchicalDirichletProcess* hdp, const char* alignment_filepath,
                               int (*kmer_to_dp_id_func) (char*), bool has_header);

#endif /* nanopore_hdp_h */
