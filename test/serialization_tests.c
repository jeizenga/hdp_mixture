//
//  serialization_tests.c
//  
//
//  Created by Jordan Eizenga on 1/15/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"
#include "hdp.h"
#include "hdp_math_utils.h"
#include "nanopore_hdp.h"
#include "sonLib.h"

void add_hdp_copy_tests(CuTest* ct, HierarchicalDirichletProcess* original_hdp,
                        HierarchicalDirichletProcess* copy_hdp) {
    
    CuAssertIntEquals_Msg(ct, "struct finalized fail, check output .hdp files",
                          (int) is_structure_finalized(original_hdp), (int) is_structure_finalized(copy_hdp));
    CuAssertIntEquals_Msg(ct, "sampling gamma fail, check output .hdp files",
                          (int) is_gamma_random(original_hdp), (int) is_gamma_random(copy_hdp));
    CuAssertIntEquals_Msg(ct, "distr finalized fail, check output .hdp files",
                          (int) is_sampling_finalized(original_hdp), (int) is_sampling_finalized(copy_hdp));
    CuAssertIntEquals_Msg(ct,  "num dir proc fail, check output .hdp files",
                          (int) get_num_dir_proc(original_hdp), (int) get_num_dir_proc(copy_hdp));
    CuAssertIntEquals_Msg(ct, "depth fail, check output .hdp files",
                          (int) get_depth(original_hdp), (int) get_depth(copy_hdp));
    CuAssertIntEquals_Msg(ct, "num data fail, check output .hdp files",
                          (int) get_num_data(original_hdp), (int) get_num_data(copy_hdp));
    CuAssertIntEquals_Msg(ct, "grid length fail, check output .hdp files",
                          (int) get_grid_length(original_hdp), (int) get_grid_length(copy_hdp));
    
    CuAssertDblEquals_Msg(ct, "mu fail, check output .hdp files",
                          get_mu(original_hdp), get_mu(copy_hdp), 0.000001);
    CuAssertDblEquals_Msg(ct, "nu fail, check output .hdp files",
                           get_nu(original_hdp), get_nu(copy_hdp), 0.0000001);
    CuAssertDblEquals_Msg(ct, "alpha fail, check output .hdp files",
                          get_alpha(original_hdp), get_alpha(copy_hdp), 0.000001);
    CuAssertDblEquals_Msg(ct, "beta fail, check output .hdp files",
                          get_beta(original_hdp), get_beta(copy_hdp), 0.000001);
    
    int64_t num_data = get_num_data(original_hdp);
    int64_t grid_length = get_grid_length(original_hdp);
    int64_t num_dps = get_num_dir_proc(original_hdp);
    
    double* original_data = get_data_copy(original_hdp);
    double* copy_data = get_data_copy(copy_hdp);
    
    int64_t* original_dp_ids = get_data_pt_dp_ids_copy(original_hdp);
    int64_t* copy_dp_ids = get_data_pt_dp_ids_copy(copy_hdp);
    
    for (int64_t i = 0; i < num_data; i++) {
        CuAssertIntEquals_Msg(ct, "data pt dp id fail, check output .hdp files",
                              (int) original_dp_ids[i], (int) copy_dp_ids[i]);
        CuAssertDblEquals_Msg(ct, "data pt fail, check output .hdp files",
                              original_data[i], copy_data[i], 0.000001);
    }
    
    free(original_data);
    free(copy_data);
    free(original_dp_ids);
    free(copy_dp_ids);
    
    int64_t original_num_dps;
    int64_t* original_num_dp_fctrs;
    int64_t original_num_gamma_params;
    double* original_gamma_params;
    double original_log_likelihood;
    double original_log_density;
    
    int64_t copy_num_dps;
    int64_t* copy_num_dp_fctrs;
    int64_t copy_num_gamma_params;
    double* copy_gamma_params;
    double copy_log_likelihood;
    double copy_log_density;
    
    take_snapshot(original_hdp, &original_num_dp_fctrs, &original_num_dps, &original_gamma_params,
                  &original_num_gamma_params, &original_log_likelihood, &original_log_density);
    
    take_snapshot(copy_hdp, &copy_num_dp_fctrs, &copy_num_dps, &copy_gamma_params,
                  &copy_num_gamma_params, &copy_log_likelihood, &copy_log_density);
    
    
    CuAssertIntEquals_Msg(ct, "num dps snapshot fail, check output .hdp files",
                          (int) original_num_dps, (int) copy_num_dps);
    CuAssertIntEquals_Msg(ct, "num gam params snapshot fail, check output .hdp files",
                          (int) original_num_gamma_params, (int) original_num_gamma_params);
    
    for (int64_t i = 0; i < original_num_dps; i++) {
        CuAssertIntEquals_Msg(ct, "num dp fctrs fail, check output .hdp files",
                              (int) original_num_dp_fctrs[i], (int) copy_num_dp_fctrs[i]);
    }
    
    for (int64_t i = 0; i < original_num_gamma_params; i++) {
        CuAssertDblEquals_Msg(ct, "gamma param fail, check output .hdp files",
                              original_gamma_params[i], copy_gamma_params[i], 0.000001);
    }
    
    CuAssertDblEquals_Msg(ct, "log likelihood fail, check output .hdp files",
                          original_log_likelihood, copy_log_likelihood, 0.000001);
    
    CuAssertDblEquals_Msg(ct, "log density fail, check output .hdp files",
                          original_log_density, copy_log_density, 0.000001);
    
    free(original_num_dp_fctrs);
    free(original_gamma_params);
    free(copy_num_dp_fctrs);
    free(copy_gamma_params);
    
    double* original_grid = get_sampling_grid_copy(original_hdp);
    double* copy_grid = get_sampling_grid_copy(copy_hdp);
    
    for (int64_t i = 0; i < grid_length; i++) {
        CuAssertDblEquals_Msg(ct, "grid fail, check output .hdp files",
                              original_grid[i], copy_grid[i], 0.000001);
    }
    
    if (is_sampling_finalized(original_hdp) && is_sampling_finalized(copy_hdp)) {
        int64_t test_grid_length = grid_length / 2;
        double* test_grid = linspace(1.1 * original_grid[0],
                                     1.1 * original_grid[grid_length - 1], test_grid_length);
        
        for (int64_t i = 0; i < test_grid_length; i++) {
            for (int64_t id = 0; id < num_dps; id++) {
                CuAssertDblEquals_Msg(ct, "dp density fail, check output .hdp files",
                                      dir_proc_density(original_hdp, test_grid[i], id),
                                      dir_proc_density(copy_hdp, test_grid[i], id),
                                      0.000001);
            }
        }
        
        free(test_grid);
    }
    
    free(original_grid);
    free(copy_grid);
    
    for (int64_t id = 0; id < num_dps; id++) {
        CuAssertIntEquals_Msg(ct, "dp size fail, check output .hdp files",
                              (int) get_dir_proc_num_factors(original_hdp, id),
                              (int) get_dir_proc_num_factors(copy_hdp, id));
        CuAssertIntEquals_Msg(ct, "dp parent fail, check output .hdp files",
                              (int) get_dir_proc_parent_id(original_hdp, id),
                              (int) get_dir_proc_parent_id(copy_hdp, id));
    }
}


void test_serialization(CuTest* ct) {
    
    FILE* data_file = fopen("/Users/Jordan/Documents/GitHub/hdp_mixture/test/data.txt","r");
    FILE* dp_id_file = fopen("/Users/Jordan/Documents/GitHub/hdp_mixture/test/dps.txt", "r");
    
    stList* data_list = stList_construct3(0, free);
    stList* dp_id_list = stList_construct3(0, free);
    
    char* data_line = stFile_getLineFromFile(data_file);
    char* dp_id_line = stFile_getLineFromFile(dp_id_file);
    
    double* datum_ptr;
    int64_t* dp_id_ptr;
    while (data_line != NULL) {
        
        datum_ptr = (double*) malloc(sizeof(double));
        dp_id_ptr = (int64_t*) malloc(sizeof(int64_t));
        
        sscanf(data_line, "%lf", datum_ptr);
        sscanf(dp_id_line, "%"SCNd64, dp_id_ptr);
        
        if (*dp_id_ptr != 4) {
            stList_append(data_list, datum_ptr);
            stList_append(dp_id_list, dp_id_ptr);
        }
        
        free(data_line);
        data_line = stFile_getLineFromFile(data_file);
        
        free(dp_id_line);
        dp_id_line = stFile_getLineFromFile(dp_id_file);
    }
    
    int64_t data_length;
    int64_t dp_ids_length;
    
    double* data = stList_toDoublePtr(data_list, &data_length);
    int64_t* dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);
    
    fclose(data_file);
    fclose(dp_id_file);
        
    int64_t num_dir_proc = 8;
    int64_t depth = 3;
    
    double mu = 0.0;
    double nu = 1.0;
    double alpha = 2.0;
    double beta = 10.0;
    
    int64_t grid_length = 250;
    double grid_start = -10.0;
    double grid_end = 10.0;
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;
    
    HierarchicalDirichletProcess* original_hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
                                                            gamma_beta, grid_start, grid_end,
                                                            grid_length, mu, nu, alpha, beta);

    
    set_dir_proc_parent(original_hdp, 1, 0);
    set_dir_proc_parent(original_hdp, 2, 0);
    set_dir_proc_parent(original_hdp, 3, 1);
    set_dir_proc_parent(original_hdp, 4, 1);
    set_dir_proc_parent(original_hdp, 5, 1);
    set_dir_proc_parent(original_hdp, 6, 2);
    set_dir_proc_parent(original_hdp, 7, 2);
    finalize_hdp_structure(original_hdp);
    
    char* filepath = "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test.hdp";
    char* copy_filepath = "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_copy.hdp";
    
    FILE* main_file;
    FILE* copy_file;
    
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    HierarchicalDirichletProcess* copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    pass_data_to_hdp(original_hdp, data, dp_ids, data_length);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    execute_gibbs_sampling(original_hdp, 10, 10, 10, false);
    finalize_distributions(original_hdp);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    destroy_hier_dir_proc(original_hdp);
    
    data = stList_toDoublePtr(data_list, &data_length);
    dp_ids = stList_toIntPtr(dp_id_list, &dp_ids_length);
    
    stList_destruct(dp_id_list);
    stList_destruct(data_list);
    
    double* gamma_params = (double*) malloc(sizeof(double) * depth);
    gamma_params[0] = 1.0; gamma_params[1] = 1.0; gamma_params[2] = 2.0;
    
    original_hdp = new_hier_dir_proc(num_dir_proc, depth, gamma_params, grid_start,
                                     grid_end, grid_length, mu, nu, alpha, beta);
    
    
    set_dir_proc_parent(original_hdp, 1, 0);
    set_dir_proc_parent(original_hdp, 2, 0);
    set_dir_proc_parent(original_hdp, 3, 1);
    set_dir_proc_parent(original_hdp, 4, 1);
    set_dir_proc_parent(original_hdp, 5, 1);
    set_dir_proc_parent(original_hdp, 6, 2);
    set_dir_proc_parent(original_hdp, 7, 2);
    finalize_hdp_structure(original_hdp);
    
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    pass_data_to_hdp(original_hdp, data, dp_ids, data_length);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    execute_gibbs_sampling(original_hdp, 10, 10, 10, false);
    finalize_distributions(original_hdp);
    main_file = fopen(filepath, "w");
    serialize_hdp(original_hdp, main_file);
    fclose(main_file);
    main_file = fopen(filepath, "r");
    copy_hdp = deserialize_hdp(main_file);
    fclose(main_file);
    copy_file = fopen(copy_filepath, "w");
    serialize_hdp(copy_hdp, copy_file);
    fclose(copy_file);
    add_hdp_copy_tests(ct, original_hdp, copy_hdp);
    destroy_hier_dir_proc(copy_hdp);
    remove(copy_filepath);
    
    destroy_hier_dir_proc(original_hdp);
    
    remove(filepath);
}

void test_nhdp_serialization(CuTest* ct) {
    
    NanoporeHDP* nhdp = flat_hdp_model("ACGT", 4, 6, 4.0, 20.0, 0.0, 100.0, 100,
                                       "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    execute_nhdp_gibbs_sampling(nhdp, 10, 0, 1, false);
    finalize_nhdp_distributions(nhdp);
    
    serialize_nhdp(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test.nhdp");
    NanoporeHDP* copy_nhdp = deserialize_nhdp("/Users/Jordan/Documents/GitHub/hdp_mixture/test/test.nhdp");
    
    
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, 65.0, "AAAAAA"),
                          get_nanopore_kmer_density(copy_nhdp, 65.0, "AAAAAA"),
                          0.000001);
    
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, 65.0, "GCACCA"),
                          get_nanopore_kmer_density(copy_nhdp, 65.0, "GCACCA"),
                          0.000001);
    
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, 65.0, "CCTTAG"),
                          get_nanopore_kmer_density(copy_nhdp, 65.0, "CCTTAG"),
                          0.000001);
    
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, 65.0, "ACTTCA"),
                          get_nanopore_kmer_density(copy_nhdp, 65.0, "ACTTCA"),
                          0.000001);
    
    CuAssertDblEquals_Msg(ct, "nhdp dp density fail",
                          get_nanopore_kmer_density(nhdp, 65.0, "GGAATC"),
                          get_nanopore_kmer_density(copy_nhdp, 65.0, "GGAATC"),
                          0.000001);
    
    remove("/Users/Jordan/Documents/GitHub/hdp_mixture/test/test.nhdp");
}

CuSuite* get_suite() {
    
    CuSuite* suite = CuSuiteNew();
    
    SUITE_ADD_TEST(suite, test_serialization);
    SUITE_ADD_TEST(suite, test_nhdp_serialization);

    return suite;
}


int main(void) {
    CuSuite* suite = get_suite();
    
    CuString* output = CuStringNew();
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    
    return 0;
}
