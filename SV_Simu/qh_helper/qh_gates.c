#include "qh_gates.h"


void applyGateSingleUnitary(Statevector *sv, int k, ComplexMatrix2 u) {
    double x_real, x_imag, y_real, y_imag;
    size_t K = 1 << k;
    double ar = u.real[0][0], br = u.real[0][1], cr = u.real[1][0], dr = u.real[1][1];
    double ai = u.imag[0][0], bi = u.imag[0][1], ci = u.imag[1][0], di = u.imag[1][1];

    for (size_t t = 0; t < sv->namp; t += (2*K))
    {
        for (size_t tt = 0; tt < K; tt++)
        {
            x_real = (ar * sv->real[t+tt] - ai * sv->imag[t+tt]) + (br * sv->real[t+tt+K] - bi * sv->imag[t+tt+K]);
            x_imag = (ar * sv->imag[t+tt] + ai * sv->real[t+tt]) + (br * sv->imag[t+tt+K] + bi * sv->real[t+tt+K]);
            y_real = (cr * sv->real[t+tt] - ci * sv->imag[t+tt]) + (dr * sv->real[t+tt+K] - di * sv->imag[t+tt+K]);
            y_imag = (cr * sv->imag[t+tt] + ci * sv->real[t+tt]) + (dr * sv->imag[t+tt+K] + di * sv->real[t+tt+K]);
            sv->real[t+tt] = x_real;
            sv->imag[t+tt] = x_imag;
            sv->real[t+tt+K] = y_real;
            sv->imag[t+tt+K] = y_imag;
        }
    }
}


void applyGateX(Statevector *sv, int k);


