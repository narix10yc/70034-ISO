#ifndef COMPLEX_MUL_H
#define COMPLEX_MUL_H

#include <stddef.h>
#include <complex.h>

typedef struct Cdouble {
    double real, imag;
} Cdouble;

void scalar_separate(double *xr, double *xi, double *yr, double *yi, double *rstr, double *rsti, size_t n);
void scalar_combined(double *x, double *y, double *rst, size_t n);
void scalar_cdouble(Cdouble *x, Cdouble *y, Cdouble *rst, size_t n);
void scalar_complexh(complex double *x, complex double *y, complex double *rst, size_t n);


void simd_separate(double *xr, double *xi, double *yr, double *yi, double *rstr, double *rsti, size_t n);
void simd_combined(double *x, double *y, double *rst, size_t n);
void simd_cdouble(Cdouble *x, Cdouble *y, Cdouble *rst, size_t n);


#endif