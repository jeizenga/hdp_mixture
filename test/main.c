#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "hdp.h"
#include "ranlib.h"

int main() {
    /*
     * This vignette walks through creating a hierarchical Dirichlet process
     * with this structure:
     *           _
     *          |H|    base normal-inverse gamma distribution (not a DP)
     *          |
     *          |_            _______
     *          |0| <--------|gamma_0|
     *        __|__
     *      _|     |_         _______
     *     |1|     |2| <-----|gamma_1|
     *  ___|___     _|_
     * |_  |_  |_  |_  |_     _______
     * |3| |4| |5| |6| |7| <-|gamma_2|
     */
	
    // the Dirichlet processes are the numbered boxes
    int num_dir_proc = 8;
    // the depth of the tree
    int depth = 3;
    // the innovation parameters at each depth (gamma and alpha_0 in Teh et al)
    double* gamma = (double*) malloc(sizeof(double) * depth);
    gamma[0] = 5.0; gamma[1] = 5.0; gamma[2] = 20.0;

    // parameters of the normal inverse gamma  base distribution
    double mu = 0.0;
    double nu = 10.0;
    double alpha = 10.0; // note: this parameter must be integer or half-integer valued
    double beta = 10.0;

    // the grid along which the HDP will record distribution samples
    int grid_length = 250;
    double* sampling_grid = NULL;// = linspace(-10.0, 10.0, grid_length);

    // initialize an HDP
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dir_proc, depth, gamma,
                                                          sampling_grid, grid_length,
                                                          mu, nu, alpha, beta);

    // establish the topology of the tree
    set_dir_proc_parent(hdp, 1, 0);
    set_dir_proc_parent(hdp, 2, 0);
    set_dir_proc_parent(hdp, 3, 1);
    set_dir_proc_parent(hdp, 4, 1);
    set_dir_proc_parent(hdp, 5, 2);
    set_dir_proc_parent(hdp, 6, 2);
    set_dir_proc_parent(hdp, 7, 2);
    // note: the result must be a perfectly balanced tree (same depth at every leaf)
    // or you will get an error here
    finalize_hdp_structure(hdp);

    //TODO: provide example data
    int data_length = 100;
    double* data = (double*) malloc(sizeof(double) * data_length); 
    int* data_pt_dps = (int*) malloc(sizeof(int) * data_length);
    for (int i = 0; i < data_length; i++) {
    	data_pt_dps[i] = 3 + (i % 5);
    	data[i] = gennor(0.0, 1.0);
    }

    pass_data_to_hdp(hdp, data, data_pt_dps, data_length);
    // note: you can also pass data before finalizing the structure
    // note: it is not necessary to observe every Dirichlet process in the data

    int num_samples = 1000;
    int burn_in = 10000;
    int thinning = 50;

    // sample from the posterior distribution of distributions
    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning);

    // calculate the mean a posteriori estimate of each distribution
    finalize_distributions(hdp);

    // query with new values
    double x = 2.1;
    int dp_id = 4;
    double density = dir_proc_density(hdp, x, dp_id);

    // reset the HDP without needing to re-initialize it and provide new data
    reset_hdp_data(hdp);

    int new_data_length = data_length;
    double* new_data = data;
    int* new_data_pt_dps = data_pt_dps;

    pass_data_to_hdp(hdp, new_data, new_data_pt_dps, new_data_length);
    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning);
    finalize_distributions(hdp);

    // query density values with the new distributions
    density = dir_proc_density(hdp, 3.4, 6);

    // free the memory
    destroy_hier_dir_proc(hdp);

    return 0;
}
