#include "two_qubit_general.h"
#include "omp.h"

void general_two_qubit_gate(Statevector *sv, int k, int l, ComplexMatrix4 u) {
    size_t K = 1 << k;
    size_t L = 1 << l;

    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];


    #pragma omp parallel \
        private(amp_pt_r, amp_pt_i, amp_r, amp_i) \
        shared(K, L, u, sv)

     #pragma omp for simd collapse(3)
    for (size_t t = 0; t < sv->namp; t += (K << 1))
    {
        for (size_t tt = 0; tt < K; tt += (L << 1))
        {
            for (size_t ttt = 0; ttt < L; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[L]   - u.imag[0][1] * amp_pt_i[L]) +
                           (u.real[0][2] * amp_pt_r[K]   - u.imag[0][2] * amp_pt_i[K]) +
                           (u.real[0][3] * amp_pt_r[L|K] - u.imag[0][3] * amp_pt_i[L|K]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[L]   - u.imag[1][1] * amp_pt_i[L]) +
                           (u.real[1][2] * amp_pt_r[K]   - u.imag[1][2] * amp_pt_i[K]) +
                           (u.real[1][3] * amp_pt_r[L|K] - u.imag[1][3] * amp_pt_i[L|K]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[L]   - u.imag[2][1] * amp_pt_i[L]) +
                           (u.real[2][2] * amp_pt_r[K]   - u.imag[2][2] * amp_pt_i[K]) +
                           (u.real[2][3] * amp_pt_r[L|K] - u.imag[2][3] * amp_pt_i[L|K]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[L]   - u.imag[3][1] * amp_pt_i[L]) +
                           (u.real[3][2] * amp_pt_r[K]   - u.imag[3][2] * amp_pt_i[K]) +
                           (u.real[3][3] * amp_pt_r[L|K] - u.imag[3][3] * amp_pt_i[L|K]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[L]   + u.imag[0][1] * amp_pt_r[L]) +
                           (u.real[0][2] * amp_pt_i[K]   + u.imag[0][2] * amp_pt_r[K]) +
                           (u.real[0][3] * amp_pt_i[L|K] + u.imag[0][3] * amp_pt_r[L|K]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[L]   + u.imag[1][1] * amp_pt_r[L]) +
                           (u.real[1][2] * amp_pt_i[K]   + u.imag[1][2] * amp_pt_r[K]) +
                           (u.real[1][3] * amp_pt_i[L|K] + u.imag[1][3] * amp_pt_r[L|K]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[L]   + u.imag[2][1] * amp_pt_r[L]) +
                           (u.real[2][2] * amp_pt_i[K]   + u.imag[2][2] * amp_pt_r[K]) +
                           (u.real[2][3] * amp_pt_i[L|K] + u.imag[2][3] * amp_pt_r[L|K]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[L]   + u.imag[3][1] * amp_pt_r[L]) +
                           (u.real[3][2] * amp_pt_i[K]   + u.imag[3][2] * amp_pt_r[K]) +
                           (u.real[3][3] * amp_pt_i[L|K] + u.imag[3][3] * amp_pt_r[L|K]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[L] = amp_r[1];
                amp_pt_r[K] = amp_r[2];
                amp_pt_r[L|K] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[L] = amp_i[1];
                amp_pt_i[K] = amp_i[2];
                amp_pt_i[L|K] = amp_i[3];
            } 
        }
    }
}


