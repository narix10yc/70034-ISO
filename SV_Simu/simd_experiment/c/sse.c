#include <smmintrin.h>
#include "simd.h"

#pragma GCC target("sse2")
void sse_add(float *x, float *y, float *rst, size_t n) {
    __m128 x_simd, y_simd, r_simd;

    for (size_t i = 0; i < n; i += 4)
    {
        x_simd = _mm_load_ps(x+i);
        y_simd = _mm_load_ps(y+i);

        r_simd = _mm_add_ps(x_simd, y_simd);

        _mm_store_ps(rst+i, r_simd);
    }
}


void sse_add_d(double *x, double *y, double *rst, size_t n) {
    __m128d x_simd, y_simd, r_simd;

    for (size_t i = 0; i < n; i += 2)
    {
        x_simd = _mm_load_pd(x+i);
        y_simd = _mm_load_pd(y+i);

        r_simd = _mm_add_pd(x_simd, y_simd);

        _mm_store_pd(rst+i, r_simd);
    }
}
