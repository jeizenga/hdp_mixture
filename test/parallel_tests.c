//
//  parallel_tests.c
//  
//
//  Created by Jordan Eizenga on 2/2/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "CuTest.h"
#include "hdp_math_utils.h"


void test_parallel_cdf(CuTest* ct) {
    
    double probs1[] = {0.2, 0.3, .45, .67, .87, .12, .322, .324, .45, .12, .56, .12, .23, .35, .39};
    double probs2[] = {0.2, 0.3, .45, .9, .67, .87, .12, .322, .324, .45, .12, .56, .12, .23, .35, .39};
    
    int64_t length1 = 15;
    int64_t length2 = 16;
    
    double* cdf = (double*) malloc(sizeof(double) * length2);
    double* cdf_direct = (double*) malloc(sizeof(double) * length2);
    
    
    double cumul = 0.0;
    for (int64_t i = 0; i < length2; i++) {
        cumul += probs1[i];
        cdf_direct[i] = cumul;
    }
    
    parallel_cdf(cdf, probs1, length1, 5);
    
    for (int64_t i = 0; i < length1; i++) {
        if (fabs(cdf[i] - cdf_direct[i]) > 0.000000001) {
            printf("round 1 failure on index %lld\n", i);
        }
        CuAssertDblEquals_Msg(ct, "cdf 1 fail\n",  cdf_direct[i], cdf[i], 0.000000001);
    }
    
    cumul = 0.0;
    for (int64_t i = 0; i < length2; i++) {
        cumul += probs2[i];
        cdf_direct[i] = cumul;
    }
    
    parallel_cdf(cdf, probs2, length2, 5);
    
    for (int64_t i = 0; i < length2; i++) {
        if (fabs(cdf[i] - cdf_direct[i]) > 0.000000001) {
            printf("round 2 failure on index %lld\n", i);
        }
        CuAssertDblEquals_Msg(ct, "cdf 2 fail\n", cdf_direct[i], cdf[i], 0.000000001);
    }
}

void test_parallel_add(CuTest* ct) {
    
    double x[] = {0.2, 0.3, .45, .67, .87, .12, .322, .324, .45, .12, .56, .12, .23, .35, .39};
    
    int64_t length = 15;
    
    double* y_direct = (double*) malloc(sizeof(double) * length);
    
    double add_factor = 4.0;
    
    for (int64_t i = 0; i < length; i++) {
        y_direct[i] = x[i] + add_factor;
    }
    
    parallel_add(add_factor, x, length);
    
    for (int64_t i = 0; i < length; i++) {
        CuAssertDblEquals_Msg(ct, "add fail\n",  y_direct[i], x[i], 0.000000001);
    }
}



void test_parallel_max(CuTest* ct) {
    
    double x[] = {0.2, 0.3, .45, .67, .87, .12, .322, .324, .45, .12, .56, .12, .23, .35, .39};
    
    int64_t length = 15;
    
    double direct_max = -1.0;
    
    for (int64_t i = 0; i < length; i++) {
        if (x[i] > direct_max) {
            direct_max = x[i];
        }
    }
    
    double par_max = parallel_max(x, length);
    
    CuAssertDblEquals_Msg(ct, "max fail\n",  par_max, direct_max, 0.000000001);
}



void test_parallel_exp(CuTest* ct) {
    
    double x[] = {0.2, 0.3, .45, .67, .87, .12, .322, .324, .45, .12, .56, .12, .23, .35, .39};
    
    int64_t length = 15;
    
    double* y_direct = (double*) malloc(sizeof(double) * length);
    
    for (int64_t i = 0; i < length; i++) {
        y_direct[i] = exp(x[i]);
    }
    
    parallel_exp(x, length);
    
    for (int64_t i = 0; i < length; i++) {
        CuAssertDblEquals_Msg(ct, "add fail\n",  y_direct[i], x[i], 0.000000001);
    }
}


CuSuite* get_suite() {
    
    CuSuite* suite = CuSuiteNew();
    
    SUITE_ADD_TEST(suite, test_parallel_cdf);
    SUITE_ADD_TEST(suite, test_parallel_add);
    SUITE_ADD_TEST(suite, test_parallel_max);
    SUITE_ADD_TEST(suite, test_parallel_exp);
    
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

