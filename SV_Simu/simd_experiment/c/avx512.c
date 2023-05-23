#include <immintrin.h>
#include "simd.h"

void avx512_add(float *x, float *y, float *rst, size_t n) {
    __m512 x_simd, y_simd, r_simd;

    for (size_t i = 0; i < n; i += 16)
    {
        x_simd = _mm512_load_ps(x+i);
        y_simd = _mm512_load_ps(y+i);

        r_simd = _mm512_add_ps(x_simd, y_simd);

        _mm512_store_ps(rst+i, r_simd);
    }
}


void avx512_add_d(double *x, double *y, double *rst, size_t n) {
    __m512d x_simd, y_simd, r_simd;

    for (size_t i = 0; i < n; i += 8)
    {
        x_simd = _mm512_load_pd(x+i);
        y_simd = _mm512_load_pd(y+i);

        r_simd = _mm512_add_pd(x_simd, y_simd);

        _mm512_store_pd(rst+i, r_simd);
    }
}


void avx512_add_d_prefetch(double *x, double *y, double *rst, size_t n) {
    __m512d x_simd, y_simd, r_simd;

    for (size_t i = 0; i < n; i += 8)
    {
        x_simd = _mm512_load_pd(x+i);
        y_simd = _mm512_load_pd(y+i);

        r_simd = _mm512_add_pd(x_simd, y_simd);

        _mm512_store_pd(rst+i, r_simd);
        _mm_prefetch(x + i + 8, _MM_HINT_T2);
        _mm_prefetch(y + i + 8, _MM_HINT_T2);

    }
}
