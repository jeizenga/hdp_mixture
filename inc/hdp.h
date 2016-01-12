#ifndef HDP_H_INCLUDED
#define HDP_H_INCLUDED

typedef struct HierarchicalDirichletProcess HierarchicalDirichletProcess;

HierarchicalDirichletProcess* new_hier_dir_proc(int num_dps, int depth, double* gamma, double sampling_grid_start,
                                                double sampling_grid_stop, int sampling_grid_length, double mu,
                                                double nu, double alpha, double beta);
HierarchicalDirichletProcess* new_hier_dir_proc_2(int num_dps, int depth, double* gamma_alpha, double* gamma_beta,
                                                  double sampling_grid_start, double sampling_grid_stop,
                                                  int sampling_grid_length, double mu, double nu, double alpha,
                                                  double beta);

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp);
void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int child_id, int parent_id);
void finalize_hdp_structure(HierarchicalDirichletProcess* hdp);
void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int* dp_id, int length);
void reset_hdp_data(HierarchicalDirichletProcess* hdp);
void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int num_samples, int burn_in, int thinning);
void execute_gibbs_sampling_with_snapshots(HierarchicalDirichletProcess* hdp,
                                           int num_samples, int burn_in, int thinning,
                                           void (*snapshot_func)(HierarchicalDirichletProcess*, void*),
                                           void* snapshot_func_args);
void finalize_distributions(HierarchicalDirichletProcess* hdp);
double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int dp_id);
void take_snapshot(HierarchicalDirichletProcess* hdp, int** num_dp_fctrs_out, int* num_dps_out,
                   double** gamma_params_out, int* num_gamma_params_out, double* log_likelihood_out);

#endif // HDP_H_INCLUDED
