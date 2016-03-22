//
//  distribution_comparison_tests.c
//  
//
//  Created by Jordan Eizenga on 1/26/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "CuTest.h"
#include "nanopore_hdp.h"
#include "hdp.h"
#include "hdp_math_utils.h"
#include "sonLib.h"

void add_distr_metric_tests(CuTest* ct, DistributionMetricMemo* memo, HierarchicalDirichletProcess* hdp) {
    int64_t num_dps = get_num_dir_proc(hdp);
    
    for (int64_t i = 0; i < num_dps; i++) {
        CuAssertDblEquals_Msg(ct, "self comparison fail\n", get_dir_proc_distance(memo, i, i),
                              0.0, 0.000000001);
        
        for (int64_t j = 0; j < i; j++) {
            double dist = get_dir_proc_distance(memo, i, j);
            
            CuAssert(ct, "nonnegativity fail\n", dist >= 0.0);
            
            CuAssertDblEquals_Msg(ct, "symmetry fail\n",  get_dir_proc_distance(memo, j, i),
                                  dist, 0.000000001);
        }
    }
}

void add_true_metric_tests(CuTest* ct, DistributionMetricMemo* memo, HierarchicalDirichletProcess* hdp) {
    int64_t num_dps = get_num_dir_proc(hdp);
    
    for (int64_t i = 0; i < num_dps - 2; i++) {
        for (int64_t j = i + 1; j < num_dps - 1; j++) {
            for (int64_t k = j + 1; k < num_dps; k++) {
                double dist_ij = get_dir_proc_distance(memo, i, j);
                double dist_jk = get_dir_proc_distance(memo, j, k);
                double dist_ik = get_dir_proc_distance(memo, i, k);
                
                bool triangle_ineq = dist_ij + dist_jk >= dist_ik - .0001;
                if (!triangle_ineq) {
                    fprintf(stderr, "failed on distances %lf + %lf < %lf\n", dist_ij, dist_jk, dist_ik);
                }
                CuAssert(ct, "triangle inequality fail\n", triangle_ineq);
                
            }
        }
    }
    
}

void test_distr_metrics(CuTest* ct) {
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
    
    stList_destruct(dp_id_list);
    stList_destruct(data_list);
    
    fclose(data_file);
    fclose(dp_id_file);
    
    int64_t num_dir_proc = 8;
    int64_t depth = 3;
    
    double mu = 0.0;
    double nu = 1.0;
    double alpha = 2.0;
    double beta = 10.0;
    
    int64_t grid_length = 500;
    double grid_start = -30.0;
    double grid_end = 30.0;
    
    double* gamma_alpha = (double*) malloc(sizeof(double) * depth);
    gamma_alpha[0] = 1.0; gamma_alpha[1] = 1.0; gamma_alpha[2] = 2.0;
    double* gamma_beta = (double*) malloc(sizeof(double) * depth);
    gamma_beta[0] = 0.2; gamma_beta[1] = 0.2; gamma_beta[2] = 0.1;
    
    HierarchicalDirichletProcess* hdp = new_hier_dir_proc_2(num_dir_proc, depth, gamma_alpha,
                                                                     gamma_beta, grid_start, grid_end,
                                                                     grid_length, mu, nu, alpha, beta);
    
    
    set_dir_proc_parent(hdp, 1, 0);
    set_dir_proc_parent(hdp, 2, 0);
    set_dir_proc_parent(hdp, 3, 1);
    set_dir_proc_parent(hdp, 4, 1);
    set_dir_proc_parent(hdp, 5, 1);
    set_dir_proc_parent(hdp, 6, 2);
    set_dir_proc_parent(hdp, 7, 2);
    finalize_hdp_structure(hdp);
    
    pass_data_to_hdp(hdp, data, dp_ids, data_length);
    
    
    execute_gibbs_sampling(hdp, 10, 10, 10, false);
    
    finalize_distributions(hdp);
    
    DistributionMetricMemo* memo = new_kl_divergence_memo(hdp);
    add_distr_metric_tests(ct, memo, hdp);
    
    memo = new_l2_distance_memo(hdp);
    add_distr_metric_tests(ct, memo, hdp);
    add_true_metric_tests(ct, memo, hdp);
    
    memo = new_shannon_jensen_distance_memo(hdp);
    add_distr_metric_tests(ct, memo, hdp);
    add_true_metric_tests(ct, memo, hdp);
    
    memo = new_hellinger_distance_memo(hdp);
    add_distr_metric_tests(ct, memo, hdp);
    add_true_metric_tests(ct, memo, hdp);
    
    destroy_hier_dir_proc(hdp);
}

void test_nhdp_distrs(CuTest* ct) {
    
    NanoporeHDP* nhdp = flat_hdp_model("ACGT", 4, 6, 4.0, 20.0, 0.0, 100.0, 100,
                                       "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_model.model");
    
    update_nhdp_from_alignment(nhdp, "/Users/Jordan/Documents/GitHub/hdp_mixture/test/test_alignment.tsv",
                               false);
    
    execute_nhdp_gibbs_sampling(nhdp, 100, 0, 1, false);
    finalize_nhdp_distributions(nhdp);
    
    NanoporeDistributionMetricMemo* memo = new_nhdp_kl_divergence_memo(nhdp);
    
    CuAssertDblEquals_Msg(ct, "kmer symmetry fail\n",  get_kmer_distr_distance(memo, "ACCCAA", "ATGATT"),
                          get_kmer_distr_distance(memo, "ATGATT", "ACCCAA"), 0.000000001);
    CuAssertDblEquals_Msg(ct, "kmer symmetry fail\n",  get_kmer_distr_distance(memo, "GCACAT", "GGGGTA"),
                          get_kmer_distr_distance(memo, "GGGGTA", "GCACAT"), 0.000000001);
    
    destroy_nanopore_hdp(nhdp);
}


CuSuite* get_suite() {
    
    CuSuite* suite = CuSuiteNew();
    
    SUITE_ADD_TEST(suite, test_distr_metrics);
    SUITE_ADD_TEST(suite, test_nhdp_distrs);
    
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

