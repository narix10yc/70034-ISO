#include "simd.h"


void scalar_add(float *x, float *y, float *rst, size_t n) {
    for (size_t i = 0; i < n; i++)
    {
        rst[i] = x[i] + y[i];
    }
}


void scalar_add_d(double *x, double *y, double *rst, size_t n) {
    for (size_t i = 0; i < n; i++)
    {
        rst[i] = x[i] + y[i];
    }
}