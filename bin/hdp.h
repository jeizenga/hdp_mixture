#ifndef HDP_H_INCLUDED
#define HDP_H_INCLUDED

typedef struct HierarchicalDirichletProcess HierarchicalDirichletProcess;

HierarchicalDirichletProcess* new_hier_dir_proc(int num_dps, int depth, double* gamma, double* sampling_grid,
                                                int grid_length, double mu, double nu, double alpha, double beta);

void destroy_hier_dir_proc(HierarchicalDirichletProcess* hdp);
void set_dir_proc_parent(HierarchicalDirichletProcess* hdp, int child_id, int parent_id);
void finalize_hdp_structure(HierarchicalDirichletProcess* hdp);
void pass_data_to_hdp(HierarchicalDirichletProcess* hdp, double* data, int* dp_id, int length);
void reset_hdp_data(HierarchicalDirichletProcess* hdp);
void execute_gibbs_sampling(HierarchicalDirichletProcess* hdp, int num_samples, int burn_in, int thinning);
void finalize_distributions(HierarchicalDirichletProcess* hdp);
double dir_proc_density(HierarchicalDirichletProcess* hdp, double x, int dp_id);

#endif // HDP_H_INCLUDED
