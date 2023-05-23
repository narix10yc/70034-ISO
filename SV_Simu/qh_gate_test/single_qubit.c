#include "single_qubit.h"
#include <immintrin.h>


void extract_indices(int n, int k, size_t *a, size_t *b) {
    size_t i, j, alpha, beta;
    for (size_t t = 0; t < (1 << (n-1)); t++)
    {
        i = (t >> k) << (k + 1);
        j = t & ((1 << k) - 1);
        alpha = i | j;
        beta = alpha ^ (1 << k);

        a[t] = alpha;
        b[t] = beta;
    }
}

void apply_on_k(int n, int k, double *amp, double a, double b, double c) {
    size_t i, j, alpha, beta;
    double x, y;
    for (size_t t = 0; t < (1 << (n-1)); t++)
    {
        i = (t >> k) << (k + 1);
        j = t & ((1 << k) - 1);
        alpha = i | j;
        beta = alpha ^ (1 << k);

        x = a * amp[alpha] + b * amp[beta];
        y = b * amp[alpha] + c * amp[beta];
        amp[alpha] = x;
        amp[beta] = y;
    }
}

void apply_on_k_prompt(int n, int k, double *amp, double a, double b, double c) {
    double x, y;
    double *amp_x, *amp_y;
    size_t K = 1 << k;
    for (size_t t = 0; t < (1 << n); t += (2*K))
    {
        for (size_t tt = 0; tt < K; tt++)
        {
            amp_x = amp + t + tt;
            amp_y = amp_x + K;
            x = a * *amp_x + b * *amp_y;
            y = b * *amp_x + c * *amp_y;
            *amp_x = x;
            *amp_y = y;
        }
    }
}

void apply_on_0(int n, double *amp, double a, double b, double c) {
    double x, y;
    for (size_t t = 0; t < (1 << n); t += 2)
    {
        x = a * amp[t] + b * amp[t+1];
        y = b * amp[t] + c * amp[t+1];
        amp[t] = x;
        amp[t+1] = y;
    }
}

// #pragma GCC target("sse4")
void apply_on_1(int n, double *amp, double a, double b, double c) {
    double x, y;
    for (size_t t = 0; t < (1 << n); t += 4)
    {
        x = a * amp[t] + b * amp[t+2];
        y = b * amp[t] + c * amp[t+2];
        amp[t] = x;
        amp[t+2] = y;

        x = a * amp[t+1] + b * amp[t+3];
        y = b * amp[t+1] + c * amp[t+3];
        amp[t+1] = x;
        amp[t+3] = y;
    }
}

void apply_on_1_manual(int n, double *amp, double a, double b, double c) {
    double *amp_x, *amp_y;
    __m128d amp_simd_x, amp_simd_y, x, y;
    __m128d a_simd, b_simd, c_simd;
    a_simd = _mm_load1_pd(&a);
    b_simd = _mm_load1_pd(&b);
    c_simd = _mm_load1_pd(&c);

    for (size_t t = 0; t < (1 << n); t += 4)
    {
        amp_x = amp + t;
        amp_y = amp + t + 2;
        amp_simd_x = _mm_load_pd(amp_x);
        amp_simd_y = _mm_load_pd(amp_y);

        x = _mm_add_pd(_mm_mul_pd(a_simd, amp_simd_x), _mm_mul_pd(b_simd, amp_simd_y));
        y = _mm_add_pd(_mm_mul_pd(b_simd, amp_simd_x), _mm_mul_pd(c_simd, amp_simd_y));
        _mm_store_pd(amp_x, x);
        _mm_store_pd(amp_y, y);
    }
}

// #pragma GCC target("avx512f")
void apply_on_2(int n, double *amp, double a, double b, double c) {
    double x, y;
    for (size_t t = 0; t < (1 << n); t += 8)
    {
        for (size_t tt = 0; tt < 4; tt++)
        {
            x = a * amp[t+tt] + b * amp[t+tt+4];
            y = b * amp[t+tt] + c * amp[t+tt+4];
            amp[t+tt] = x;
            amp[t+tt+4] = y;
        }
    }
}

// #pragma GCC target("avx")
void apply_on_2_manual(int n, double *amp, double a, double b, double c) {
    double *amp_x, *amp_y;
    __m256d amp_simd_x, amp_simd_y, x, y;
    __m256d a_simd, b_simd, c_simd;
    a_simd = _mm256_set1_pd(a);
    b_simd = _mm256_set1_pd(b);
    c_simd = _mm256_set1_pd(c);

    for (size_t t = 0; t < (1 << n); t += 8)
    {
        amp_x = amp + t;
        amp_y = amp + t + 4;
        amp_simd_x = _mm256_load_pd(amp_x);
        amp_simd_y = _mm256_load_pd(amp_y);

        x = _mm256_add_pd(_mm256_mul_pd(a_simd, amp_simd_x), _mm256_mul_pd(b_simd, amp_simd_y));
        y = _mm256_add_pd(_mm256_mul_pd(b_simd, amp_simd_x), _mm256_mul_pd(c_simd, amp_simd_y));
        _mm256_store_pd(amp_x, x);
        _mm256_store_pd(amp_y, y);
    }
}

void apply_on_3(int n, double *amp, double a, double b, double c) {
    double x, y;
    for (size_t t = 0; t < (1 << n); t += 16)
    {
        for (size_t tt = 0; tt < 8; tt++)
        {
            x = a * amp[t+tt] + b * amp[t+tt+8];
            y = b * amp[t+tt] + c * amp[t+tt+8];
            amp[t+tt] = x;
            amp[t+tt+8] = y;
        }
    }
}

// #pragma GCC target("avx512f")
void apply_on_3_manual(int n, double *amp, double a, double b, double c) {
    double *amp_x, *amp_y;
    __m512d amp_simd_x, amp_simd_y, x, y;
    __m512d a_simd, b_simd, c_simd;
    a_simd = _mm512_set1_pd(a);
    b_simd = _mm512_set1_pd(b);
    c_simd = _mm512_set1_pd(c);

    for (size_t t = 0; t < (1 << n); t += 16)
    {
        amp_x = amp + t;
        amp_y = amp_x + 8;
        amp_simd_x = _mm512_load_pd(amp_x);
        amp_simd_y = _mm512_load_pd(amp_y);

        x = _mm512_add_pd(_mm512_mul_pd(a_simd, amp_simd_x), _mm512_mul_pd(b_simd, amp_simd_y));
        y = _mm512_add_pd(_mm512_mul_pd(b_simd, amp_simd_x), _mm512_mul_pd(c_simd, amp_simd_y));
        _mm512_store_pd(amp_x, x);
        _mm512_store_pd(amp_y, y);
    }
}

void apply_on_4(int n, double *amp, double a, double b, double c) {
    double x, y;
    for (size_t t = 0; t < (1 << n); t += 32)
    {
        for (size_t tt = 0; tt < 16; tt++)
        {
            x = a * amp[t+tt] + b * amp[t+tt+16];
            y = b * amp[t+tt] + c * amp[t+tt+16];
            amp[t+tt] = x;
            amp[t+tt+16] = y;
        }
    }
}

// #pragma GCC target("avx512f")
void apply_on_4_manual(int n, double *amp, double a, double b, double c) {
    double *amp_x, *amp_y;
    __m512d amp_simd_x, amp_simd_y, x, y;
    __m512d a_simd, b_simd, c_simd;
    a_simd = _mm512_set1_pd(a);
    b_simd = _mm512_set1_pd(b);
    c_simd = _mm512_set1_pd(c);

    for (size_t t = 0; t < (1 << n); t += 32)
    {
        for (size_t tt = 0; tt < 16; tt += 8)
        {
            amp_x = amp + t + tt;
            amp_y = amp_x + 16;
            amp_simd_x = _mm512_load_pd(amp_x);
            amp_simd_y = _mm512_load_pd(amp_y);

            x = _mm512_add_pd(_mm512_mul_pd(a_simd, amp_simd_x), _mm512_mul_pd(b_simd, amp_simd_y));
            y = _mm512_add_pd(_mm512_mul_pd(b_simd, amp_simd_x), _mm512_mul_pd(c_simd, amp_simd_y));
            _mm512_store_pd(amp_x, x);
            _mm512_store_pd(amp_y, y);
        }
    }
}