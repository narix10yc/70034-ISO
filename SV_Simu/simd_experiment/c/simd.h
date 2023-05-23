#ifndef SIMD_H
#define SIMD_H

#include <stddef.h>

void no_sse_scalar_add(float *x, float *y, float *rst, size_t n);

void scalar_add(float *x, float *y, float *rst, size_t n);

void sse_add(float *x, float *y, float *rst, size_t n);

void avx_add(float *x, float *y, float *rst, size_t n);

void avx512_add(float *x, float *y, float *rst, size_t n);

void no_sse_scalar_add_d(double *x, double *y, double *rst, size_t n);

void scalar_add_d(double *x, double *y, double *rst, size_t n);

void sse_add_d(double *x, double *y, double *rst, size_t n);

void avx_add_d(double *x, double *y, double *rst, size_t n);

void avx512_add_d(double *x, double *y, double *rst, size_t n);

void avx512_add_d_prefetch(double *x, double *y, double *rst, size_t n);

#endif