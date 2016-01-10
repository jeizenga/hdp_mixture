#ifndef HDP_MATH_UTILS_H_INCLUDED
#define HDP_MATH_UTILS_H_INCLUDED

#include <stdbool.h>

typedef struct SumOfLogsMemo SumOfLogsMemo;

SumOfLogsMemo* new_log_sum_memo();
void destroy_log_sum_memo(SumOfLogsMemo* memo);

// returns log(Gamma(n / 2)) in amortized constant time with low risk of overflow
double log_gamma_half(int n, SumOfLogsMemo* sum_of_logs_memo);
double sum_of_logs(SumOfLogsMemo* memo, int n);

// quick-select on array copy (does not alter original array)
double median(double* arr, int length);

// returns the index of the first element of arr greater or equal to x, assuming arr is sorted
// returns final index if x is greater than all elements of arr
int bisect_left(double x, double* arr, int length);

double* spline_knot_slopes(double* x, double* y, int length);
double spline_interp(double query_x, double* x, double* y, double* slope, int length);
double grid_spline_interp(double query_x, double* x, double* y, double* slope, int length);

double double_max(double a, double b);

double* linspace(double start, double stop, int length);

double rand_uniform(double a);
double rand_beta(double a, double b);
bool rand_bernoulli(double p);

// explained further in Jordan's math notebook section "Cached variables for improved performance"
double log_posterior_conditional_term(double nu_post, double two_alpha_post, double beta_post, SumOfLogsMemo* memo);

#endif // HDP_MATH_UTILS_H_INCLUDED
