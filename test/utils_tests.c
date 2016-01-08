#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "hdp.h"
#include "hdp_math_utils.h"


bool close(double a, double b) {
    return (a - b > .00000001) && (b - a > .00000001);
}

int main()
{

    int arr_len = 10;
    double* arr = (double*) malloc(sizeof(double) * arr_len);
    for (int i = 0; i < arr_len; i++) {
        arr[i] = (double) i;
    }

    // bisect tests

    for (int i = 0; i < arr_len - 1; i++) {
        if (bisect_left((double) i - 0.25, arr, arr_len) != i) {
            printf("bisect test 1 failed on %d", i);
            exit(EXIT_FAILURE);
        }
        if (bisect_left((double) i, arr, arr_len) != i) {
            printf("bisect test 2 failed on %d", i);
            exit(EXIT_FAILURE);
        }
        if (bisect_left((double) i + 0.25, arr, arr_len) != i + 1) {
            printf("bisect test 3 failed on %d", i);
            exit(EXIT_FAILURE);
        }
    }

    if (bisect_left(6.0, arr, arr_len) != arr_len - 1) {
        printf("bisect test 4 failed");
        exit(EXIT_FAILURE);
    }

    // median tests

    double med;

    med = median(arr, arr_len);
    if (med != (double) (arr_len / 2)) {
        printf("median test 1 failed");
        exit(EXIT_FAILURE);
    }

    arr[0] = (double) arr_len;
    med = median(arr, arr_len);
    if (med != (double) (arr_len / 2 + 1)) {
        printf("median test 2 failed");
        exit(EXIT_FAILURE);
    }

    arr[arr_len - 1] = (double) arr_len;
    med = median(arr, arr_len);
    if (med != (double) (arr_len / 2 + 1)) {
        printf("median test 3 failed");
        exit(EXIT_FAILURE);
    }

    arr[1] = 2.0 * (double) arr_len;
    med = median(arr, arr_len);
    if (med != (double) (arr_len / 2 + 2)) {
        printf("median test 3 failed");
        exit(EXIT_FAILURE);
    }

    free(arr);

    int length = 5;
    double* x = linspace(0.0, 1.0, length);

    if (!close(x[0], 0.0)) {
        printf("linspace test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(x[1], 0.25)) {
        printf("linspace test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(x[2], 0.5)) {
        printf("linspace test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(x[3], 0.75)) {
        printf("linspace test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(x[4], 1.0)) {
        printf("linspace test failed");
        exit(EXIT_FAILURE);
    }

    double* y = (double*) malloc(sizeof(double) * length);
    y[0] = 1.0;
    y[1] = 4.0;
    y[2] = 0.0;
    y[3] = 3.0;
    y[4] = -5.5;

    double* slopes = spline_knot_slopes(x, y, length);

    if (!close(slopes[0], 22.3214285714)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(slopes[1], -8.64285714286)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(slopes[2], 0.25)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(slopes[3], -4.35714285714)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
    if (!close(slopes[4], -48.8214285714)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
	
    double* query_x = linspace(-0.1, 1.1, 6);
    double y_interp;

    for (int i = 0; i < length; i++) {
        if (!close(spline_interp(x[i], x, y, slopes, length), y[i])) {
            printf("spline slope basic test failed");
            exit(EXIT_FAILURE);
        }
    }

    y_interp = spline_interp(query_x[0], x, y, slopes, length);
    if (!close(y_interp, -1.2321428571428572)) {
        printf("spline slope test failed");
        exit(EXIT_FAILURE);
    }
    y_interp = spline_interp(query_x[1], x, y, slopes, length);
    if (!close(y_interp, 3.6718480000000007)) {
        printf("spline interp test failed");
        exit(EXIT_FAILURE);
    }
    y_interp = spline_interp(query_x[2], x, y, slopes, length);
    if (!close(y_interp, 1.6130811428571403)) {
        printf("spline interp test failed");
        exit(EXIT_FAILURE);
    }
    y_interp = spline_interp(query_x[3], x, y, slopes, length);
    if (!close(y_interp, 1.548665142857147)) {
        printf("spline interp test failed");
        exit(EXIT_FAILURE);
    }
    y_interp = spline_interp(query_x[4], x, y, slopes, length);
    if (!close(y_interp, 0.68427999999999289)) {
        printf("spline interp test failed");
        exit(EXIT_FAILURE);
    }
    y_interp = spline_interp(query_x[5], x, y, slopes, length);
    if (!close(y_interp, -10.382142857142862)) {
        printf("spline interp test failed");
        exit(EXIT_FAILURE);
    }

    free(query_x);
    free(x);
    free(y);
    free(slopes);

    SumOfLogsMemo* memo = new_log_sum_memo();
    //for (int i = 1; i < 10; i++) {
    //    double manual_sum = 0.0;
    //    for (int j = 1; j <= i; j++) {
    //        manual_sum += log((double) j);
    //    }
    //    if (!close(sum_of_logs(memo, i), manual_sum)) {
    //        printf("log sum test failed");
    //        exit(EXIT_FAILURE);
    //    }
    //}

    for (int i = 2; i < 20; i += 3) {
        if (!close(log_gamma_half(i, memo), lgamma(.5 * (double) i))) {
            printf("gamma function test failed");
            exit(EXIT_FAILURE);
        }
    }

    double alpha = 10.0;
    double beta = 10.0;
    double nu = 11.0;
    double lpct = log_posterior_conditional_term(nu, 2.0 * alpha, beta, memo);
    if (!close(lpct, -11.422971086258172)) {
        printf("log posterior conditional function test 1 failed");
        exit(EXIT_FAILURE);
    }

    alpha = 5.0;
    beta = 11.0;
    nu = 9.0;
    lpct = log_posterior_conditional_term(nu, 2.0 * alpha, beta, memo);
    if (!close(lpct, -9.9100348223120)) {
        printf("log posterior conditional function test 2 failed");
        exit(EXIT_FAILURE);
    }

    alpha = 31.1;
    beta = 8.5;
    nu = 100.78;
    lpct = log_posterior_conditional_term(nu, 2.0 * alpha, beta, memo);
    if (!close(lpct / 100000.0, 0.0000613765)) {
        printf("log posterior conditional function test 3 failed");
        exit(EXIT_FAILURE);
    }

    destroy_log_sum_memo(memo);

    return 0;
}
