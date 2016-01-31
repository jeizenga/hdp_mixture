//
//  base_params_tests.c
//  
//
//  Created by Jordan Eizenga on 1/30/16.
//
//

#include <stdio.h>
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
