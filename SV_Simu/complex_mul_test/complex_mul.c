#include "complex_mul.h"
#include <immintrin.h>

void scalar_separate(double *xr, double *xi, double *yr, double *yi, double *rstr, double *rsti, size_t n) {
    for (size_t i = 0; i < n; i ++)
    {
        rstr[i] = (xr[i] * yr[i]) - (xi[i] * yi[i]);
        rsti[i] = (xr[i] * yi[i]) + (xi[i] * yr[i]);
    }
}

void scalar_combined(double *x, double *y, double *rst, size_t n) {
    for (size_t i = 0; i < n; i ++)
    {
        rst[2*i] = (x[2*i] * y[2*i]) - (x[2*i+1] * y[2*i+1]);
        rst[2*i+1] = (x[2*i] * y[2*i+1]) + (x[2*i+1] * y[2*i]);
    }
}

void scalar_cdouble(Cdouble *x, Cdouble *y, Cdouble *rst, size_t n) {
    for (size_t i = 0; i < n; i ++)
    {
        rst[i].real = (x[i].real * y[i].real) - (x[i].imag * y[i].imag);
        rst[i].imag = (x[i].real * y[i].imag) + (x[i].imag * y[i].real);
    }
}

void scalar_complexh(complex double *x, complex double *y, complex double *rst, size_t n) {
    for (size_t i = 0; i < n; i++)
    {
        rst[i] = x[i] * y[i];
    }
    
}


void simd_separate(double *xr, double *xi, double *yr, double *yi, double *rstr, double *rsti, size_t n) {
    __m512d simd_xr, simd_xi, simd_yr, simd_yi, simd_rstr, simd_rsti;
    for (size_t t = 0; t < n; t += 8)
    {
        simd_xr = _mm512_load_pd(xr+t);
        simd_xi = _mm512_load_pd(xi+t);
        simd_yr = _mm512_load_pd(yr+t);
        simd_yi = _mm512_load_pd(yi+t);

        simd_rstr = _mm512_fmsub_pd(simd_xr, simd_yr, _mm512_mul_pd(simd_xi, simd_yi));
        simd_rsti = _mm512_fmadd_pd(simd_xr, simd_yi, _mm512_mul_pd(simd_xi, simd_yr));

        _mm512_store_pd(rstr+t, simd_rstr);
        _mm512_store_pd(rsti+t, simd_rsti);
    }
}


void simd_combined(double *x, double *y, double *rst, size_t n) {
    __m512d simd_x, simd_y, zmm0, zmm1, zmm2, zmm3, zmm_real, zmm_imag;
    for (size_t t = 0; t < 2*n; t += 8)
    {
        simd_x = _mm512_load_pd(x+t); // xr xi
        simd_y = _mm512_load_pd(y+t); // yr yi

        zmm0 = _mm512_mul_pd(simd_x, simd_y); // xr*yr xi*yi
        simd_y = _mm512_shuffle_pd(simd_y, simd_y, 0b01010101); // yi yr

        zmm1 = _mm512_mul_pd(simd_x, simd_y); // xr*yi xi*yr

        zmm2 = _mm512_shuffle_pd(zmm0, zmm1, 0b00000000); // xi*yr xi*yi
        zmm3 = _mm512_shuffle_pd(zmm0, zmm1, 0b11111111); // xr*yi xr*yr

        zmm_imag = _mm512_add_pd(zmm2, zmm3); // imag part
        zmm_real = _mm512_add_pd(zmm3, zmm2); // real part

        _mm512_mask_store_pd(rst+t, 0b10101010, zmm_real);
        _mm512_mask_store_pd(rst+t, 0b01010101, zmm_real);

    }
    

}


void simd_cdouble(Cdouble *x, Cdouble *y, Cdouble *rst, size_t n) {

}