#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "hdp.h"
#include "sonLib.h"
#include "hdp_math_utils.h"

void load_data(const char* data_filepath, const char* dp_id_filepath,
               double** data_out, int64_t** dp_ids_out, int64_t* length_out) {
    
    FILE* data_file = fopen(data_filepath, "r");
    FILE* dp_id_file = fopen(dp_id_filepath, "r");
    
    stList* data_list = stList_construct3(0, *free);
    stList* dp_id_list = stList_construct3(0, *free);
    
    char* line = stFile_getLineFromFile(data_file);
    
    double* datum_ptr;
    while (line != NULL) {

        datum_ptr = (double*) malloc(sizeof(double));
        
        sscanf(line, "%lf", datum_ptr);
        
        stList_append(data_list, datum_ptr);
        
        free(line);
        line = stFile_getLineFromFile(data_file);
    }
    
    
    line = stFile_getLineFromFile(dp_id_file);
    int64_t* dp_id_ptr;
    while (line != NULL) {
        dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));
        
        sscanf(line, "%"SCNd64, dp_id_ptr);
        
        stList_append(dp_id_list, dp_id_ptr);
        
        free(line);
        line = stFile_getLineFromFile(dp_id_file);
    }
    
    int64_t data_length;
    int64_t dp_ids_length;
    
    double* data = stList_toDoublePtr(data_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);
    
    if (data_length != dp_ids_length) {
        fprintf(stderr, "Data and DP ID files have different lengths\n");
        exit(EXIT_FAILURE);
    }
    
    *data_out = data;
    *dp_ids_out = dp_ids;
    *length_out = data_length;
    
    stList_destruct(dp_id_list);
    stList_destruct(data_list);
    
    fclose(data_file);
    fclose(dp_id_file);
}

void write_double_data(FILE* f, double* data, int64_t length) {
    if (length <= 0) {
        return;
    }
    fprintf(f, "%lf", data[0]);
    for (int i = 1; i < length; i++) {
        fprintf(f, "\n%lf", data[i]);
    }
}

void write_int_data(FILE* f, int64_t* data, int64_t length) {
    if (length <= 0) {
        return;
    }
    fprintf(f, "%"PRId64, data[0]);
    for (int64_t i = 1; i < length; i++) {
        fprintf(f, "\n%"PRId64, data[i]);
    }
}

void output_distrs_to_disk(HierarchicalDirichletProcess* hdp, double* grid, int64_t grid_length, int64_t num_dps) {
    FILE* x = fopen("x_vals.txt", "w");
    write_double_data(x, grid, grid_length);
    fclose(x);
    
    char file_name[30];
    
    double* pdf = (double*) malloc(sizeof(double) * grid_length);
    for (int64_t i = 0; i < num_dps; i++) {
        for (int64_t j = 0; j < grid_length; j++) {
            pdf[j] = dir_proc_density(hdp, grid[j], i);
        }
        sprintf(file_name, "dp_%"PRId64"_distr.txt", i);
        FILE* out = fopen(file_name, "w");
        write_double_data(out, pdf, grid_length);
        fclose(out);
    }
    free(pdf);
}

typedef struct SnapshotArgs {
    const char* num_dp_factors_filepath;
    const char* gamma_params_filepath;
    const char* log_likelihood_filepath;
} SnapshotArgs;

SnapshotArgs make_snapshot_args(const char* num_dp_factors_filepath,
                                const char* gamma_params_filepath,
                                const char* log_likelihood_filepath) {
    SnapshotArgs args;
    args.num_dp_factors_filepath = num_dp_factors_filepath;
    args.gamma_params_filepath = gamma_params_filepath;
    args.log_likelihood_filepath = log_likelihood_filepath;
    return args;
}

void record_snapshots_to_files(HierarchicalDirichletProcess* hdp, void* snapshot_args_void_ptr) {
    SnapshotArgs* snapshot_args_ptr = (SnapshotArgs*) snapshot_args_void_ptr;
    SnapshotArgs snapshot_args = *snapshot_args_ptr;
    
    FILE* num_dp_factors_file = fopen(snapshot_args.num_dp_factors_filepath, "a");
    FILE* gamma_params_file = fopen(snapshot_args.gamma_params_filepath, "a");
    FILE* log_likelihood_file = fopen(snapshot_args.log_likelihood_filepath, "a");
    
    int64_t* num_dp_factors;
    int64_t num_dps;
    
    double* gamma_params;
    int64_t num_gamma_params;
    
    double log_likelihood;
    
    take_snapshot(hdp, &num_dp_factors, &num_dps, &gamma_params, &num_gamma_params, &log_likelihood);
    
    fprintf(log_likelihood_file, "%lf\n", log_likelihood);
    
    for (int64_t i = 0; i < num_gamma_params - 1; i++) {
        fprintf(gamma_params_file, "%lf\t", gamma_params[i]);
    }
    fprintf(gamma_params_file, "%lf\n", gamma_params[num_gamma_params - 1]);
    
    for (int64_t i = 0; i < num_dps - 1; i++) {
        fprintf(num_dp_factors_file, "%"PRId64"\t", num_dp_factors[i]);
    }
    fprintf(num_dp_factors_file, "%"PRId64"\n", num_dp_factors[num_dps - 1]);
    
    fclose(num_dp_factors_file);
    fclose(gamma_params_file);
    fclose(log_likelihood_file);
}



int main(int argc, char* argv[]) {
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
    
    if (argc > 2) {
        fprintf(stderr, "Too many command line arguments\n");
    }
    else if (argc == 2) {
        srand(atoi(argv[1]));
    }
	
    // the Dirichlet processes are the numbered boxes
    int64_t num_dir_proc = 8;
    // the depth of the tree
    int64_t depth = 3;

    // parameters of the normal inverse gamma  base distribution
    double mu = 0.0;
    double nu = 1.0;
    double alpha = 2.0; // note: this parameter must be integer or half-integer valued
    double beta = 10.0;

    // the grid along which the HDP will record distribution samples
    int64_t grid_length = 250;
    double grid_start = -10.0;
    double grid_end = 10.0;

    // initialize an HDP
    fprintf(stderr, "Initializing HDP...\n");
    
    // choose whether to pre-define concentration parameters (gamma and alpha_0 in Teh et al)
    // or sample them from a Gamma distribution
    
    // pre-define the concentration parameters at each depth
    double* gamma = (double*) malloc(sizeof(double) * depth);
    gamma[0] = 10.0; gamma[1] = 10.0; gamma[2] = 10.0;
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dir_proc, depth, gamma,
                                                          grid_start, grid_end, grid_length,
                                                          mu, nu, alpha, beta);
    
    
    // parameters for distributions of concentration parameters at each depth
//    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
//    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
//    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
//    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;
//    HierarchicalDirichletProcess* hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
//                                                            gamma_beta, grid_start, grid_end,
//                                                            grid_length, mu, nu, alpha, beta);
    
    // establish the topology of the tree
    fprintf(stderr, "Establishing HDP tree topology...\n");
    set_dir_proc_parent(hdp, 1, 0);
    set_dir_proc_parent(hdp, 2, 0);
    set_dir_proc_parent(hdp, 3, 1);
    set_dir_proc_parent(hdp, 4, 1);
    set_dir_proc_parent(hdp, 5, 1);
    set_dir_proc_parent(hdp, 6, 2);
    set_dir_proc_parent(hdp, 7, 2);
    
//    // the Dirichlet processes are the numbered boxes
//    int64_t num_dir_proc = 21;
//    // the depth of the tree
//    int64_t depth = 2;
//    
//    // parameters of the normal inverse gamma  base distribution
//    double mu = 0.0;
//    double nu = 1.0;
//    double alpha = 1.5; // note: this parameter must be integer or half-integer valued
//    double beta = 10.0;
//    
//    // the grid along which the HDP will record distribution samples
//    int64_t grid_length = 250;
//    double grid_start = -10.0;
//    double grid_end = 10.0;
//    
//    // initialize an HDP
//    fprintf(stderr, "Initializing HDP...\n");
//    
//    // choose whether to pre-define concentration parameters (gamma and alpha_0 in Teh et al)
//    // or sample them from a Gamma distribution
//    
//    //    // pre-define the concentration parameters at each depth
//        double* gamma = (double*) malloc(sizeof(double) * depth);
//        gamma[0] = 10.0; gamma[1] = 5.0;
//        HierarchicalDirichletProcess* hdp = new_hier_dir_proc(num_dir_proc, depth, gamma,
//                                                              grid_start, grid_end, grid_length,
//                                                              mu, nu, alpha, beta);
    
    
    // parameters for distributions of concentration parameters at each depth
//    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
//    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
//    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
//    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;
//    HierarchicalDirichletProcess* hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
//                                                            gamma_beta, grid_start, grid_end,
//                                                            grid_length, mu, nu, alpha, beta);
//    
//    // establish the topology of the tree
//    fprintf(stderr, "Establishing HDP tree topology...\n");
//    for (int64_t dp_id = 0; dp_id < 20; dp_id++) {
//        set_dir_proc_parent(hdp, dp_id, 20);
//    }
    
    // note: the result must be a perfectly balanced tree (same depth at every leaf)
    // or you will get an error here
    finalize_hdp_structure(hdp);
    
    fprintf(stderr, "Loading data from disk...\n");
    double* data;
    int64_t* data_pt_dps;
    int64_t data_length;
    load_data("/Users/Jordan/Documents/GitHub/hdp_mixture/test/data.txt",
              "/Users/Jordan/Documents/GitHub/hdp_mixture/test/dps.txt",
              &data, &data_pt_dps, &data_length);
    
    int64_t new_data_length = data_length / 2;
    double* new_data = (double*) malloc(sizeof(double) * new_data_length);
    int64_t* new_data_pt_dps = (int64_t*) malloc(sizeof(int64_t) * new_data_length);
    
    for (int64_t i = 0; i < new_data_length; i++) {
        new_data[i] = data[i] + 1.0;
        new_data_pt_dps[i] = data_pt_dps[i];
    }
    
    fprintf(stderr, "Giving HDP data...\n");
    pass_data_to_hdp(hdp, data, data_pt_dps, data_length);
    // note: you can also pass data before finalizing the structure
    // note: it is not necessary to observe every Dirichlet process in the data
    
    int64_t num_samples = 2500;
    int64_t burn_in = 2000000;
    int64_t thinning = 500;

    // choose whether to Gibbs sample only the distributions or to also supply a
    // snapshot function that samples at the beginning of each Gibbs sweep
    
//    // sample with a snapshot function
//    fprintf(stderr, "Making snapshot args...\n");
//    SnapshotArgs filepaths = make_snapshot_args("num_dp_factors.txt",
//                                                "gamma_params.txt",
//                                                "log_likelihood.txt");
//    
//    fprintf(stderr, "Executing Gibbs sampling and recording snapshots...\n");
//    execute_gibbs_sampling_with_snapshots(hdp, num_samples, burn_in, thinning,
//                                          &record_snapshots_to_files, (void*) &filepaths, true);
    
    // sample without a snapshot function
    fprintf(stderr, "Executing Gibbs sampling without snapshots...\n");
    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning, true);

    // calculate the mean a posteriori estimate of each distribution
    fprintf(stderr, "Computing mean a posteriori distributions...\n");
    finalize_distributions(hdp);

    // query with new values
    double x = 2.1;
    int64_t dp_id = 4;
    double density = dir_proc_density(hdp, x, dp_id);

    int64_t x_len = 200;
    double x_start = -10.0;
    double x_end = 10.0;
    double* x_vals = (double*) malloc(sizeof(double) * x_len);
    for (int64_t i = 0; i < x_len; i++) {
        x_vals[i] = x_start + ((double) i) * (x_end - x_start) / ((double) (x_len - 1));
    }

    output_distrs_to_disk(hdp, x_vals, x_len, num_dir_proc);

//    // reset the HDP without needing to re-initialize it and provide new data
//    fprintf(stderr, "Reseting HDP data...\n");
//    reset_hdp_data(hdp);
//    
//    fprintf(stderr, "Giving HDP new data...\n");
//    pass_data_to_hdp(hdp, new_data, new_data_pt_dps, new_data_length);
//    fprintf(stderr, "Executing Gibbs sampling without snapshots...\n");
//    execute_gibbs_sampling(hdp, num_samples, burn_in, thinning, true);
//    fprintf(stderr, "Computing mean a posteriori distributions...\n");
//    finalize_distributions(hdp);
//
//    // query density values with the new distributions
//    double new_density = dir_proc_density(hdp, 3.4, 6);

    // free the memory
    fprintf(stderr, "Destroying HDP...\n");
    destroy_hier_dir_proc(hdp);
    
    printf("Completed!\n");

    return 0;
}
