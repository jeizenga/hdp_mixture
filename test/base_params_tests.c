//
//  base_params_tests.c
//  
//
//  Created by Jordan Eizenga on 1/30/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "CuTest.h"
#include "hdp_math_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338
#endif

double norm_gamma_density(double mu, double tau, double mu_0, double nu, double alpha, double beta) {
    return pow(beta, alpha) / tgamma(alpha) * pow(tau, alpha - 1.0) * exp(- beta * tau)
           * sqrt(nu * tau / (2.0 * M_PI)) * exp(-(nu * tau / 2.0) * pow(mu - mu_0, 2.0));
}

double norm_gamma_joint_log_likelihood(double* mus, double* taus, int64_t length,
                                       double mu_0, double nu, double alpha, double beta) {
    double log_likelihood = 0.0;
    for (int64_t i = 0; i < length; i++) {
        log_likelihood += log(norm_gamma_density(mus[i], taus[i], mu_0, nu, alpha, beta));
    }
    
    return log_likelihood;
}

double t_predictive_log_likelihood(double x, double mu_0, double nu, double alpha, double beta) {
    return lgamma(alpha + .5) - lgamma(alpha) + .5 * log(nu / (2.0 * M_PI *(nu + 1.0) * beta))
           - (alpha + .5) * log(1.0 + nu / (2.0 *(nu + 1.0) * beta) * pow(x - mu_0, 2.0));
}

double t_predictive_joint_log_likelihood(double* xs, int64_t length,
                                         double mu_0, double nu, double alpha, double beta) {
    double log_likelihood = 0.0;
    for (int64_t i = 0; i < length; i++) {
        log_likelihood += t_predictive_log_likelihood(xs[i], mu_0, nu, alpha, beta);
    }
    
    return log_likelihood;
}

void test_lineq_solve(CuTest* ct) {
    //printf("lineq\n");
    static double A[] = {1.0, 2.0, 3.0, 2.0, 3.0, 1.0, 0.0, 1.0, 1.0};
    static double b[] = {14.0, 11.0, 5.0};
    static double x[] = {1.0, 2.0, 3.0};
    
    double* solution = (double*) malloc(sizeof(double) * 3);
    
    lineq_solve(A, b, solution, 3);
    
    for (int64_t i = 0; i < 3; i++) {
        CuAssertDblEquals_Msg(ct, "lin eq solve fail.\n",
                              x[i], solution[i], 0.000001);
    }
    
}

void test_inverse(CuTest* ct) {
    //printf("inverse\n");
    static double A[] = {2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0};
    static double B[] = {1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0};
    static double C[] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    
    double* A_inv = matrix_inverse(A, 3);
    double* B_inv = matrix_inverse(B, 3);
    double* C_inv = matrix_inverse(C, 3);
    
    
    for (int64_t i = 0; i < 3; i++) {
        for (int64_t j = 0; j < 3; j++) {
            if (i == j) {
                CuAssertDblEquals_Msg(ct, "matrix inverse fail.\n",
                                      A_inv[i * 3 + j], 0.5, 0.000001);
            }
            else {
                CuAssertDblEquals_Msg(ct, "matrix inverse fail.\n",
                                     A_inv[i * 3 + j], 0.0, 0.000001);
            }
        }
    }
    
    for (int64_t i = 0; i < 3; i++) {
        for (int64_t j = 0; j < 3; j++) {
            if (i <= j) {
                CuAssertDblEquals_Msg(ct, "matrix inverse fail.\n",
                                      B_inv[i * 3 + j], 1.0, 0.000001);
            }
            else {
                CuAssertDblEquals_Msg(ct, "matrix inverse fail.\n",
                                      B_inv[i * 3 + j], 0.0, 0.000001);
            }
        }
    }
    
    for (int64_t i = 0; i < 3; i++) {
        for (int64_t j = 0; j < 3; j++) {
            CuAssertDblEquals_Msg(ct, "matrix inverse fail.\n",
                                  C_inv[i * 3 + j], C[i * 3 + j], 0.000001);
        }
    }
    
    free(A_inv);
    free(B_inv);
    free(C_inv);
    
}

void test_matrix_mult(CuTest* ct) {
    //printf("mult\n");
    static double A[] = {2.0, 1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0, 2.0};
    static double b[] = {14.0, 11.0, 5.0};
    static double B[] = {14.0, 11.0, 5.0, 6.0, 4.0, 1.0};
    double* c = matrix_mult(A, b, 3, 3, 1);
    double* C = matrix_mult(A, B, 3, 3, 2);
    
    
    for (int64_t i = 0; i < 3; i++) {
        CuAssertDblEquals_Msg(ct, "matrix mult fail.\n",
                              c[i], 2.0 * b[i] + b[(i + 1) % 3], 0.000001);
    }
    
    for (int64_t i = 0; i < 3; i++) {
        for (int64_t j = 0; i < 2; i++) {
            CuAssertDblEquals_Msg(ct, "matrix mult fail.\n",
                                  C[i * 2 + j], 2.0 * B[i * 2 + j] + B[((i + 1) % 3) * 2 + j], 0.000001);
        }
    }
    
    free(c);
    free(C);
}

void test_max_pred_params(CuTest* ct) {
    //printf("max pred\n");
    //static double x[] = {-20.1, 2.8, -11.7, -39.3, -0.4};
    static double x[] = {0.0, 0.0, 0.0, 0.0, 1.0};
    int64_t length = 5;
    
    double mu_0, nu, alpha, beta;
    //max_pred_normal_inverse_gamma_params(x, length, &mu_0, &nu, &alpha, &beta);
    max_pred_normal_inverse_gamma_params_2(x, length, 0.0001, &mu_0, &nu, &alpha, &beta, 0.0000000001);
    
    double max_likelihood = t_predictive_joint_log_likelihood(x, length, mu_0, nu, alpha, beta);
    
    double factor = 1.414;
    int64_t min_power = -4;
    int64_t max_power = 5;
    
    for (int64_t i = min_power; i < max_power; i++) {
        for (int64_t j = min_power; j < max_power; j++) {
            for (int64_t k = min_power; k < max_power; k++) {
                for (int64_t l = min_power; l < max_power; l++) {
                    double candidate_likelihood = t_predictive_joint_log_likelihood(x, length,
                                                                                    pow(factor, i) * mu_0,
                                                                                    pow(factor, j) * nu,
                                                                                    pow(factor, k) * alpha,
                                                                                    pow(factor, l) * beta);
                    
                    CuAssert(ct, "higher likelihood exists\n",
                             candidate_likelihood <= max_likelihood + .0000001);
                }
            }
        }
    }
}

void test_mle_params(CuTest* ct) {
    static double mus[] = {-20.1, 2.8, -11.7, -39.3, -0.4};
    static double taus[] = {0.01, 0.005, 0.0023, 0.013, 0.008};
    int64_t length = 5;
    
    double mu_0, nu, alpha, beta;
    
    mle_normal_inverse_gamma_params(mus, taus, length, &mu_0, &nu, &alpha, &beta);
    
    double max_likelihood = norm_gamma_joint_log_likelihood(mus, taus, length, mu_0, nu, alpha, beta);
    //printf("####### max: %lf\n", max_likelihood);
    
    for (int64_t i = -2; i < 3; i++) {
        for (int64_t j = -2; j < 3; j++) {
            for (int64_t k = -2; k < 3; k++) {
                for (int64_t l = -2; l < 3; l++) {
                    double candidate_likelihood = norm_gamma_joint_log_likelihood(mus, taus, length,
                                                                                  pow(2.0, i) * mu_0,
                                                                                  pow(2.0, j) * nu,
                                                                                  pow(2.0, k) * alpha,
                                                                                  pow(2.0, l) * beta);
                    
                    //printf("%lf\n", candidate_likelihood);
                    CuAssert(ct, "higher likelihood exists\n",
                             candidate_likelihood <= max_likelihood + .0000001);
                }
            }
        }
    }
}

CuSuite* get_suite() {
    
    CuSuite* suite = CuSuiteNew();
    
    SUITE_ADD_TEST(suite, test_mle_params);
    SUITE_ADD_TEST(suite, test_lineq_solve);
    SUITE_ADD_TEST(suite, test_inverse);
    SUITE_ADD_TEST(suite, test_matrix_mult);
    SUITE_ADD_TEST(suite, test_max_pred_params);
    
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
