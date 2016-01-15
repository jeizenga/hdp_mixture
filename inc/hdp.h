#ifndef HDP_H_INCLUDED
#define HDP_H_INCLUDED

typedef struct HierarchicalDirichletProcess HierarchicalDirichletProcess;

HierarchicalDirichletProcess* new_hier_dir_proc(int64_t num_dps, int64_t depth, double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop, int64_t sampling_grid_length, double mu,
                                                double nu, double alpha, double beta);
HierarchicalDirichletProcess* new_hier_dir_proc_2(int64_t num_dps, int64_t depth, double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start, double sampling_grid_stop,
                                                  int64_t sampling_grid_length, double mu, double nu, double alpha,
                                                  double beta);

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp);

void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int64_t child_id, int64_t parent_id);

void finalize_hdp_structure(HierarchicalDirichletProcess* hdp);

void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int64_t* dp_id, int64_t length);

void reset_hdp_data(HierarchicalDirichletProcess* hdp);

void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int64_t num_samples, int64_t burn_in, int64_t thinning);

void execute_gibbs_sampling_with_snapshots(HierarchicalDirichletProcess* hdp,
                                           int64_t num_samples, int64_t burn_in, int64_t thinning,
                                           void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                           void* snapshot_func_args);

void finalize_distributions(HierarchicalDirichletProcess* hdp);

double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int64_t dp_id);

void take_snapshot(HierarchicalDirichletProcess* hdp, int64_t** num_dp_fctrs_out, int64_t* num_dps_out,
                   double** gamma_params_out, int64_t* num_gamma_params_out, double* log_likelihood_out);

#endif // HDP_H_INCLUDED
