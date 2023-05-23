#include <immintrin.h>
#include <omp.h>
#include <stdio.h>

#include "single_qubit.h"


// void extract_indices(int n, int k, size_t *a, size_t *b) {
//     size_t i, j, alpha, beta;
//     for (size_t t = 0; t < (1 << (n-1)); t++)
//     {
//         i = (t >> k) << (k + 1);
//         j = t & ((1 << k) - 1);
//         alpha = i | j;
//         beta = alpha ^ (1 << k);

//         a[t] = alpha;
//         b[t] = beta;
//     }
// }

void apply_on_k(int n, int k, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    size_t i, j, alpha, beta; // Get indices
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << (n-1)); t++)
    {
        i = (t >> k) << (k + 1);
        j = t & ((1 << k) - 1);
        alpha = i | j;
        beta = alpha ^ (1 << k);

        x_real = (a.real * amp.real[alpha] - a.imag * amp.imag[alpha]) + (b.real * amp.real[beta] - b.imag * amp.imag[beta]);
        x_imag = (a.real * amp.imag[alpha] + a.imag * amp.real[alpha]) + (b.real * amp.imag[beta] + b.imag * amp.real[beta]);
        y_real = (c.real * amp.real[alpha] - c.imag * amp.imag[alpha]) + (d.real * amp.real[beta] - d.imag * amp.imag[beta]);
        y_imag = (c.real * amp.imag[alpha] + c.imag * amp.real[alpha]) + (d.real * amp.imag[beta] + d.imag * amp.real[beta]);
        amp.real[alpha] = x_real;
        amp.imag[alpha] = x_imag;
        amp.real[beta] = y_real;
        amp.imag[beta] = y_imag;
    }
}

// void apply_on_k_prompt(int n, int k, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
//     double x_real, x_imag, y_real, y_imag;
//     double *amp_x_real, *amp_x_imag, *amp_y_real, *amp_y_imag;
//     size_t K = 1 << k;
//     for (size_t t = 0; t < (1 << n); t += (2*K))
//     {
//         for (size_t tt = 0; tt < K; tt++)
//         {
//             amp_x_real = amp.real + t + tt;
//             amp_x_imag = amp.imag + t + tt;
//             amp_y_real = amp_x_real + K;
//             amp_y_imag = amp_x_imag + K;
//             x_real = (a.real * *amp_x_real - a.imag * *amp_x_imag) + (b.real * *amp_y_real - b.imag * *amp_y_imag);
//             x_imag = (a.real * *amp_x_imag + a.imag * *amp_x_real) + (b.real * *amp_y_imag + b.imag * *amp_y_real);
//             y_real = (c.real * *amp_x_real - c.imag * *amp_x_imag) + (d.real * *amp_y_real - d.imag * *amp_y_imag);
//             y_imag = (c.real * *amp_x_imag + c.imag * *amp_x_real) + (d.real * *amp_y_imag + d.imag * *amp_y_real);

//             *amp_x_real = x_real;
//             *amp_y_real = y_real;
//             *amp_x_imag = x_imag;
//             *amp_y_imag = y_imag;
//         }
//     }
// }


void apply_on_k_prompt(int n, int k, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    size_t K = 1 << k;        
    for (size_t t = 0; t < (1 << n); t += (2*K))
    {
        for (size_t tt = 0; tt < K; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+K] - b.imag * amp.imag[t+tt+K]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+K] + b.imag * amp.real[t+tt+K]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+K] - d.imag * amp.imag[t+tt+K]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+K] + d.imag * amp.real[t+tt+K]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+K] = y_real;
            amp.imag[t+tt+K] = y_imag;
            
        }
    }

}

void apply_on_k_prompt_simd(int n, int k, Camp amp, ComplexMatrix2 u) {
    double x_real, x_imag, y_real, y_imag;
    size_t K = 1 << k;
    double ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1];
    double ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1];
    #pragma omp parallel \
        // shared (amp, ar, br, cr, dr, ai, bi, ci, di, K) \
        // private (x_real, x_imag, y_real, y_imag)
    {
        if (k == 0)
        {
            #pragma omp for
            for (size_t t = 0; t < (1 << n); t += 2)
            {
                x_real = (ar * amp.real[t] - ai * amp.imag[t]) + (br * amp.real[t+1] - bi * amp.imag[t+1]);
                x_imag = (ar * amp.imag[t] + ai * amp.real[t]) + (br * amp.imag[t+1] + bi * amp.real[t+1]);
                y_real = (cr * amp.real[t] - ci * amp.imag[t]) + (dr * amp.real[t+1] - di * amp.imag[t+1]);
                y_imag = (cr * amp.imag[t] + ci * amp.real[t]) + (dr * amp.imag[t+1] + di * amp.real[t+1]);
                amp.real[t] = x_real;
                amp.imag[t] = x_imag;
                amp.real[t+1] = y_real;
                amp.imag[t+1] = y_imag;
            }
        }
        else 
        {
            #pragma omp for simd collapse(2)
            for (size_t t = 0; t < (1 << n); t += (2*K))
            {
                // #pragma omp simd
                for (size_t tt = 0; tt < K; tt++)
                {
                    x_real = (ar * amp.real[t+tt] - ai * amp.imag[t+tt]) + (br * amp.real[t+tt+K] - bi * amp.imag[t+tt+K]);
                    x_imag = (ar * amp.imag[t+tt] + ai * amp.real[t+tt]) + (br * amp.imag[t+tt+K] + bi * amp.real[t+tt+K]);
                    y_real = (cr * amp.real[t+tt] - ci * amp.imag[t+tt]) + (dr * amp.real[t+tt+K] - di * amp.imag[t+tt+K]);
                    y_imag = (cr * amp.imag[t+tt] + ci * amp.real[t+tt]) + (dr * amp.imag[t+tt+K] + di * amp.real[t+tt+K]);
                    amp.real[t+tt] = x_real;
                    amp.imag[t+tt] = x_imag;
                    amp.real[t+tt+K] = y_real;
                    amp.imag[t+tt+K] = y_imag;
                }
            }
        }
    }
}

void apply_on_k_prompt_matrix(int n, int k, Camp amp, ComplexMatrix2 u) {
    double x_real, x_imag, y_real, y_imag;
    size_t K = 1 << k;
    double ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1];
    double ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1];
    #pragma omp parallel \
        // shared (ar, br, cr, dr, ai, bi, ci, di, K) \
        // private (x_real, x_imag, y_real, y_imag)
    {
        #pragma omp for
        for (size_t t = 0; t < (1 << n); t += (2*K))
        {
            // #pragma omp for
            for (size_t tt = 0; tt < K; tt++)
            {
                x_real = (ar * amp.real[t+tt] - ai * amp.imag[t+tt]) + (br * amp.real[t+tt+K] - bi * amp.imag[t+tt+K]);
                x_imag = (ar * amp.imag[t+tt] + ai * amp.real[t+tt]) + (br * amp.imag[t+tt+K] + bi * amp.real[t+tt+K]);
                y_real = (cr * amp.real[t+tt] - ci * amp.imag[t+tt]) + (dr * amp.real[t+tt+K] - di * amp.imag[t+tt+K]);
                y_imag = (cr * amp.imag[t+tt] + ci * amp.real[t+tt]) + (dr * amp.imag[t+tt+K] + di * amp.real[t+tt+K]);
                amp.real[t+tt] = x_real;
                amp.imag[t+tt] = x_imag;
                amp.real[t+tt+K] = y_real;
                amp.imag[t+tt+K] = y_imag;
            }
        }
    }
}

void apply_on_k_prompt_multi_threading(int n, int k, Camp amp, ComplexMatrix2 u) {
    double x_real, x_imag, y_real, y_imag;
    size_t K = 1 << k;
    double ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1];
    double ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1];
    #pragma omp parallel \
        shared (ar, br, cr, dr, ai, bi, ci, di, K) \
        private (x_real, x_imag, y_real, y_imag)
    {
        if (k == 0)
        {
            #pragma omp for
            for (size_t t = 0; t < (1 << n); t += 2)
            {
                x_real = (ar * amp.real[t] - ai * amp.imag[t]) + (br * amp.real[t+1] - bi * amp.imag[t+1]);
                x_imag = (ar * amp.imag[t] + ai * amp.real[t]) + (br * amp.imag[t+1] + bi * amp.real[t+1]);
                y_real = (cr * amp.real[t] - ci * amp.imag[t]) + (dr * amp.real[t+1] - di * amp.imag[t+1]);
                y_imag = (cr * amp.imag[t] + ci * amp.real[t]) + (dr * amp.imag[t+1] + di * amp.real[t+1]);
                amp.real[t] = x_real;
                amp.imag[t] = x_imag;
                amp.real[t+1] = y_real;
                amp.imag[t+1] = y_imag;
            }
        }
        else 
        {
            #pragma omp for simd collapse(2)
            for (size_t t = 0; t < (1 << n); t += (2*K))
            {
                // #pragma omp simd
                for (size_t tt = 0; tt < K; tt++)
                {
                    x_real = (ar * amp.real[t+tt] - ai * amp.imag[t+tt]) + (br * amp.real[t+tt+K] - bi * amp.imag[t+tt+K]);
                    x_imag = (ar * amp.imag[t+tt] + ai * amp.real[t+tt]) + (br * amp.imag[t+tt+K] + bi * amp.real[t+tt+K]);
                    y_real = (cr * amp.real[t+tt] - ci * amp.imag[t+tt]) + (dr * amp.real[t+tt+K] - di * amp.imag[t+tt+K]);
                    y_imag = (cr * amp.imag[t+tt] + ci * amp.real[t+tt]) + (dr * amp.imag[t+tt+K] + di * amp.real[t+tt+K]);
                    amp.real[t+tt] = x_real;
                    amp.imag[t+tt] = x_imag;
                    amp.real[t+tt+K] = y_real;
                    amp.imag[t+tt+K] = y_imag;
                }
            }
        }
    }
}


void explicit_on_0(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 2)
    {
        x_real = (a.real * amp.real[t] - a.imag * amp.imag[t]) + (b.real * amp.real[t+1] - b.imag * amp.imag[t+1]);
        x_imag = (a.real * amp.imag[t] + a.imag * amp.real[t]) + (b.real * amp.imag[t+1] + b.imag * amp.real[t+1]);
        y_real = (c.real * amp.real[t] - c.imag * amp.imag[t]) + (d.real * amp.real[t+1] - d.imag * amp.imag[t+1]);
        y_imag = (c.real * amp.imag[t] + c.imag * amp.real[t]) + (d.real * amp.imag[t+1] + d.imag * amp.real[t+1]);
        amp.real[t] = x_real;
        amp.imag[t] = x_imag;
        amp.real[t+1] = y_real;
        amp.imag[t+1] = y_imag;
    }
}

// // #pragma GCC target("sse4")
void explicit_on_1(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 4)
    {
        x_real = (a.real * amp.real[t] - a.imag * amp.imag[t]) + (b.real * amp.real[t+2] - b.imag * amp.imag[t+2]);
        x_imag = (a.real * amp.imag[t] + a.imag * amp.real[t]) + (b.real * amp.imag[t+2] + b.imag * amp.real[t+2]);
        y_real = (c.real * amp.real[t] - c.imag * amp.imag[t]) + (d.real * amp.real[t+2] - d.imag * amp.imag[t+2]);
        y_imag = (c.real * amp.imag[t] + c.imag * amp.real[t]) + (d.real * amp.imag[t+2] + d.imag * amp.real[t+2]);
        amp.real[t] = x_real;
        amp.imag[t] = x_imag;
        amp.real[t+2] = y_real;
        amp.imag[t+2] = y_imag;

        x_real = (a.real * amp.real[t+1] - a.imag * amp.imag[t+1]) + (b.real * amp.real[t+3] - b.imag * amp.imag[t+3]);
        x_imag = (a.real * amp.imag[t+1] + a.imag * amp.real[t+1]) + (b.real * amp.imag[t+3] + b.imag * amp.real[t+3]);
        y_real = (c.real * amp.real[t+1] - c.imag * amp.imag[t+1]) + (d.real * amp.real[t+3] - d.imag * amp.imag[t+3]);
        y_imag = (c.real * amp.imag[t+1] + c.imag * amp.real[t+1]) + (d.real * amp.imag[t+3] + d.imag * amp.real[t+3]);
        amp.real[t+1] = x_real;
        amp.imag[t+1] = x_imag;
        amp.real[t+3] = y_real;
        amp.imag[t+3] = y_imag;
    }
}


void explicit_on_2(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 8)
    {
        for (size_t tt = 0; tt < 4; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+4] - b.imag * amp.imag[t+tt+4]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+4] + b.imag * amp.real[t+tt+4]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+4] - d.imag * amp.imag[t+tt+4]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+4] + d.imag * amp.real[t+tt+4]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+4] = y_real;
            amp.imag[t+tt+4] = y_imag;
        }
    }
}


void explicit_on_3(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 16)
    {
        for (size_t tt = 0; tt < 8; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+8] - b.imag * amp.imag[t+tt+8]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+8] + b.imag * amp.real[t+tt+8]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+8] - d.imag * amp.imag[t+tt+8]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+8] + d.imag * amp.real[t+tt+8]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+8] = y_real;
            amp.imag[t+tt+8] = y_imag;
        }
    }
}

void explicit_on_4(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 32)
    {
        for (size_t tt = 0; tt < 16; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+16] - b.imag * amp.imag[t+tt+16]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+16] + b.imag * amp.real[t+tt+16]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+16] - d.imag * amp.imag[t+tt+16]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+16] + d.imag * amp.real[t+tt+16]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+16] = y_real;
            amp.imag[t+tt+16] = y_imag;
        }
    }
}

void explicit_on_5(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 64)
    {
        for (size_t tt = 0; tt < 32; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+32] - b.imag * amp.imag[t+tt+32]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+32] + b.imag * amp.real[t+tt+32]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+32] - d.imag * amp.imag[t+tt+32]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+32] + d.imag * amp.real[t+tt+32]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+32] = y_real;
            amp.imag[t+tt+32] = y_imag;
        }
    }
}

void explicit_on_6(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 128)
    {
        for (size_t tt = 0; tt < 64; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+64] - b.imag * amp.imag[t+tt+64]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+64] + b.imag * amp.real[t+tt+64]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+64] - d.imag * amp.imag[t+tt+64]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+64] + d.imag * amp.real[t+tt+64]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+64] = y_real;
            amp.imag[t+tt+64] = y_imag;
        }
    }
}

void explicit_on_7(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 256)
    {
        for (size_t tt = 0; tt < 128; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+128] - b.imag * amp.imag[t+tt+128]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+128] + b.imag * amp.real[t+tt+128]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+128] - d.imag * amp.imag[t+tt+128]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+128] + d.imag * amp.real[t+tt+128]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+128] = y_real;
            amp.imag[t+tt+128] = y_imag;
        }
    }
}

void explicit_on_8(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 512)
    {
        for (size_t tt = 0; tt < 256; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+256] - b.imag * amp.imag[t+tt+256]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+256] + b.imag * amp.real[t+tt+256]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+256] - d.imag * amp.imag[t+tt+256]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+256] + d.imag * amp.real[t+tt+256]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+256] = y_real;
            amp.imag[t+tt+256] = y_imag;
        }
    }
}

void explicit_on_9(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 1024)
    {
        for (size_t tt = 0; tt < 512; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+512] - b.imag * amp.imag[t+tt+512]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+512] + b.imag * amp.real[t+tt+512]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+512] - d.imag * amp.imag[t+tt+512]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+512] + d.imag * amp.real[t+tt+512]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+512] = y_real;
            amp.imag[t+tt+512] = y_imag;
        }
    }
}

void explicit_on_10(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d) {
    double x_real, x_imag, y_real, y_imag;
    for (size_t t = 0; t < (1 << n); t += 2048)
    {
        for (size_t tt = 0; tt < 1024; tt++)
        {
            x_real = (a.real * amp.real[t+tt] - a.imag * amp.imag[t+tt]) + (b.real * amp.real[t+tt+1024] - b.imag * amp.imag[t+tt+1024]);
            x_imag = (a.real * amp.imag[t+tt] + a.imag * amp.real[t+tt]) + (b.real * amp.imag[t+tt+1024] + b.imag * amp.real[t+tt+1024]);
            y_real = (c.real * amp.real[t+tt] - c.imag * amp.imag[t+tt]) + (d.real * amp.real[t+tt+1024] - d.imag * amp.imag[t+tt+1024]);
            y_imag = (c.real * amp.imag[t+tt] + c.imag * amp.real[t+tt]) + (d.real * amp.imag[t+tt+1024] + d.imag * amp.real[t+tt+1024]);
            amp.real[t+tt] = x_real;
            amp.imag[t+tt] = x_imag;
            amp.real[t+tt+1024] = y_real;
            amp.imag[t+tt+1024] = y_imag;
        }
    }
}


// void apply_on_1_manual(int n, double *amp, double a, double b, double c) {
//     double *amp_x, *amp_y;
//     __m128d amp_simd_x, amp_simd_y, x, y;
//     __m128d a_simd, b_simd, c_simd;
//     a_simd = _mm_load1_pd(&a);
//     b_simd = _mm_load1_pd(&b);
//     c_simd = _mm_load1_pd(&c);

//     for (size_t t = 0; t < (1 << n); t += 4)
//     {
//         amp_x = amp + t;
//         amp_y = amp + t + 2;
//         amp_simd_x = _mm_load_pd(amp_x);
//         amp_simd_y = _mm_load_pd(amp_y);

//         x = _mm_add_pd(_mm_mul_pd(a_simd, amp_simd_x), _mm_mul_pd(b_simd, amp_simd_y));
//         y = _mm_add_pd(_mm_mul_pd(b_simd, amp_simd_x), _mm_mul_pd(c_simd, amp_simd_y));
//         _mm_store_pd(amp_x, x);
//         _mm_store_pd(amp_y, y);
//     }
// }



// // #pragma GCC target("avx")
// void apply_on_2_manual(int n, double *amp, double a, double b, double c) {
//     double *amp_x, *amp_y;
//     __m256d amp_simd_x, amp_simd_y, x, y;
//     __m256d a_simd, b_simd, c_simd;
//     a_simd = _mm256_set1_pd(a);
//     b_simd = _mm256_set1_pd(b);
//     c_simd = _mm256_set1_pd(c);

//     for (size_t t = 0; t < (1 << n); t += 8)
//     {
//         amp_x = amp + t;
//         amp_y = amp + t + 4;
//         amp_simd_x = _mm256_load_pd(amp_x);
//         amp_simd_y = _mm256_load_pd(amp_y);

//         x = _mm256_add_pd(_mm256_mul_pd(a_simd, amp_simd_x), _mm256_mul_pd(b_simd, amp_simd_y));
//         y = _mm256_add_pd(_mm256_mul_pd(b_simd, amp_simd_x), _mm256_mul_pd(c_simd, amp_simd_y));
//         _mm256_store_pd(amp_x, x);
//         _mm256_store_pd(amp_y, y);
//     }
// }



// // #pragma GCC target("avx512f")
// void apply_on_3_manual(int n, double *amp, double a, double b, double c) {
//     double *amp_x, *amp_y;
//     __m512d amp_simd_x, amp_simd_y, x, y;
//     __m512d a_simd, b_simd, c_simd;
//     a_simd = _mm512_set1_pd(a);
//     b_simd = _mm512_set1_pd(b);
//     c_simd = _mm512_set1_pd(c);

//     for (size_t t = 0; t < (1 << n); t += 16)
//     {
//         amp_x = amp + t;
//         amp_y = amp_x + 8;
//         amp_simd_x = _mm512_load_pd(amp_x);
//         amp_simd_y = _mm512_load_pd(amp_y);

//         x = _mm512_add_pd(_mm512_mul_pd(a_simd, amp_simd_x), _mm512_mul_pd(b_simd, amp_simd_y));
//         y = _mm512_add_pd(_mm512_mul_pd(b_simd, amp_simd_x), _mm512_mul_pd(c_simd, amp_simd_y));
//         _mm512_store_pd(amp_x, x);
//         _mm512_store_pd(amp_y, y);
//     }
// }



// // #pragma GCC target("avx512f")
// void apply_on_4_manual(int n, double *amp, double a, double b, double c) {
//     double *amp_x, *amp_y;
//     __m512d amp_simd_x, amp_simd_y, x, y;
//     __m512d a_simd, b_simd, c_simd;
//     a_simd = _mm512_set1_pd(a);
//     b_simd = _mm512_set1_pd(b);
//     c_simd = _mm512_set1_pd(c);

//     for (size_t t = 0; t < (1 << n); t += 32)
//     {
//         for (size_t tt = 0; tt < 16; tt += 8)
//         {
//             amp_x = amp + t + tt;
//             amp_y = amp_x + 16;
//             amp_simd_x = _mm512_load_pd(amp_x);
//             amp_simd_y = _mm512_load_pd(amp_y);

//             x = _mm512_add_pd(_mm512_mul_pd(a_simd, amp_simd_x), _mm512_mul_pd(b_simd, amp_simd_y));
//             y = _mm512_add_pd(_mm512_mul_pd(b_simd, amp_simd_x), _mm512_mul_pd(c_simd, amp_simd_y));
//             _mm512_store_pd(amp_x, x);
//             _mm512_store_pd(amp_y, y);
//         }
//     }
// }









// QuEST
void statevec_unitaryLocal(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    double stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    double *stateVecReal = qureg.real;
    double *stateVecImag = qureg.imag;

// // # ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u,numTasks) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo)
// // # endif
//     {
// // # ifdef _OPENMP
# pragma omp for schedule (static) 
// // # endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            // store current state vector values in temp variables
            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateRealLo = stateVecReal[indexLo];
            stateImagLo = stateVecImag[indexLo];


            // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
            stateVecReal[indexUp] = u.real[0][0]*stateRealUp - u.imag[0][0]*stateImagUp 
                + u.real[0][1]*stateRealLo - u.imag[0][1]*stateImagLo;
            stateVecImag[indexUp] = u.real[0][0]*stateImagUp + u.imag[0][0]*stateRealUp 
                + u.real[0][1]*stateImagLo + u.imag[0][1]*stateRealLo;

            // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
            stateVecReal[indexLo] = u.real[1][0]*stateRealUp  - u.imag[1][0]*stateImagUp 
                + u.real[1][1]*stateRealLo  -  u.imag[1][1]*stateImagLo;
            stateVecImag[indexLo] = u.real[1][0]*stateImagUp + u.imag[1][0]*stateRealUp 
                + u.real[1][1]*stateImagLo + u.imag[1][1]*stateRealLo;

        // } 
    }
} 
