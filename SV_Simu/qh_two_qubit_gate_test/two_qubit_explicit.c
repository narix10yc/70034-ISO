#include "two_qubit_explicit.h"


void two_qubit_gate_1_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (2 << 1))
    {
        for (size_t tt = 0; tt < 2; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[2]   - u.imag[0][2] * amp_pt_i[2]) +
                           (u.real[0][3] * amp_pt_r[3] - u.imag[0][3] * amp_pt_i[3]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[2]   - u.imag[1][2] * amp_pt_i[2]) +
                           (u.real[1][3] * amp_pt_r[3] - u.imag[1][3] * amp_pt_i[3]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[2]   - u.imag[2][2] * amp_pt_i[2]) +
                           (u.real[2][3] * amp_pt_r[3] - u.imag[2][3] * amp_pt_i[3]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[2]   - u.imag[3][2] * amp_pt_i[2]) +
                           (u.real[3][3] * amp_pt_r[3] - u.imag[3][3] * amp_pt_i[3]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[2]   + u.imag[0][2] * amp_pt_r[2]) +
                           (u.real[0][3] * amp_pt_i[3] + u.imag[0][3] * amp_pt_r[3]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[2]   + u.imag[1][2] * amp_pt_r[2]) +
                           (u.real[1][3] * amp_pt_i[3] + u.imag[1][3] * amp_pt_r[3]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[2]   + u.imag[2][2] * amp_pt_r[2]) +
                           (u.real[2][3] * amp_pt_i[3] + u.imag[2][3] * amp_pt_r[3]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[2]   + u.imag[3][2] * amp_pt_r[2]) +
                           (u.real[3][3] * amp_pt_i[3] + u.imag[3][3] * amp_pt_r[3]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[2] = amp_r[2];
                amp_pt_r[3] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[2] = amp_i[2];
                amp_pt_i[3] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_2_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (4 << 1))
    {
        for (size_t tt = 0; tt < 4; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[4]   - u.imag[0][2] * amp_pt_i[4]) +
                           (u.real[0][3] * amp_pt_r[5] - u.imag[0][3] * amp_pt_i[5]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[4]   - u.imag[1][2] * amp_pt_i[4]) +
                           (u.real[1][3] * amp_pt_r[5] - u.imag[1][3] * amp_pt_i[5]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[4]   - u.imag[2][2] * amp_pt_i[4]) +
                           (u.real[2][3] * amp_pt_r[5] - u.imag[2][3] * amp_pt_i[5]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[4]   - u.imag[3][2] * amp_pt_i[4]) +
                           (u.real[3][3] * amp_pt_r[5] - u.imag[3][3] * amp_pt_i[5]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[4]   + u.imag[0][2] * amp_pt_r[4]) +
                           (u.real[0][3] * amp_pt_i[5] + u.imag[0][3] * amp_pt_r[5]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[4]   + u.imag[1][2] * amp_pt_r[4]) +
                           (u.real[1][3] * amp_pt_i[5] + u.imag[1][3] * amp_pt_r[5]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[4]   + u.imag[2][2] * amp_pt_r[4]) +
                           (u.real[2][3] * amp_pt_i[5] + u.imag[2][3] * amp_pt_r[5]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[4]   + u.imag[3][2] * amp_pt_r[4]) +
                           (u.real[3][3] * amp_pt_i[5] + u.imag[3][3] * amp_pt_r[5]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[4] = amp_r[2];
                amp_pt_r[5] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[4] = amp_i[2];
                amp_pt_i[5] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_3_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (8 << 1))
    {
        for (size_t tt = 0; tt < 8; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[8]   - u.imag[0][2] * amp_pt_i[8]) +
                           (u.real[0][3] * amp_pt_r[9] - u.imag[0][3] * amp_pt_i[9]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[8]   - u.imag[1][2] * amp_pt_i[8]) +
                           (u.real[1][3] * amp_pt_r[9] - u.imag[1][3] * amp_pt_i[9]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[8]   - u.imag[2][2] * amp_pt_i[8]) +
                           (u.real[2][3] * amp_pt_r[9] - u.imag[2][3] * amp_pt_i[9]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[8]   - u.imag[3][2] * amp_pt_i[8]) +
                           (u.real[3][3] * amp_pt_r[9] - u.imag[3][3] * amp_pt_i[9]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[8]   + u.imag[0][2] * amp_pt_r[8]) +
                           (u.real[0][3] * amp_pt_i[9] + u.imag[0][3] * amp_pt_r[9]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[8]   + u.imag[1][2] * amp_pt_r[8]) +
                           (u.real[1][3] * amp_pt_i[9] + u.imag[1][3] * amp_pt_r[9]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[8]   + u.imag[2][2] * amp_pt_r[8]) +
                           (u.real[2][3] * amp_pt_i[9] + u.imag[2][3] * amp_pt_r[9]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[8]   + u.imag[3][2] * amp_pt_r[8]) +
                           (u.real[3][3] * amp_pt_i[9] + u.imag[3][3] * amp_pt_r[9]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[8] = amp_r[2];
                amp_pt_r[9] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[8] = amp_i[2];
                amp_pt_i[9] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_4_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (16 << 1))
    {
        for (size_t tt = 0; tt < 16; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[16]   - u.imag[0][2] * amp_pt_i[16]) +
                           (u.real[0][3] * amp_pt_r[17] - u.imag[0][3] * amp_pt_i[17]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[16]   - u.imag[1][2] * amp_pt_i[16]) +
                           (u.real[1][3] * amp_pt_r[17] - u.imag[1][3] * amp_pt_i[17]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[16]   - u.imag[2][2] * amp_pt_i[16]) +
                           (u.real[2][3] * amp_pt_r[17] - u.imag[2][3] * amp_pt_i[17]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[16]   - u.imag[3][2] * amp_pt_i[16]) +
                           (u.real[3][3] * amp_pt_r[17] - u.imag[3][3] * amp_pt_i[17]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[16]   + u.imag[0][2] * amp_pt_r[16]) +
                           (u.real[0][3] * amp_pt_i[17] + u.imag[0][3] * amp_pt_r[17]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[16]   + u.imag[1][2] * amp_pt_r[16]) +
                           (u.real[1][3] * amp_pt_i[17] + u.imag[1][3] * amp_pt_r[17]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[16]   + u.imag[2][2] * amp_pt_r[16]) +
                           (u.real[2][3] * amp_pt_i[17] + u.imag[2][3] * amp_pt_r[17]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[16]   + u.imag[3][2] * amp_pt_r[16]) +
                           (u.real[3][3] * amp_pt_i[17] + u.imag[3][3] * amp_pt_r[17]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[16] = amp_r[2];
                amp_pt_r[17] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[16] = amp_i[2];
                amp_pt_i[17] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_5_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (32 << 1))
    {
        for (size_t tt = 0; tt < 32; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[32]   - u.imag[0][2] * amp_pt_i[32]) +
                           (u.real[0][3] * amp_pt_r[33] - u.imag[0][3] * amp_pt_i[33]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[32]   - u.imag[1][2] * amp_pt_i[32]) +
                           (u.real[1][3] * amp_pt_r[33] - u.imag[1][3] * amp_pt_i[33]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[32]   - u.imag[2][2] * amp_pt_i[32]) +
                           (u.real[2][3] * amp_pt_r[33] - u.imag[2][3] * amp_pt_i[33]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[32]   - u.imag[3][2] * amp_pt_i[32]) +
                           (u.real[3][3] * amp_pt_r[33] - u.imag[3][3] * amp_pt_i[33]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[32]   + u.imag[0][2] * amp_pt_r[32]) +
                           (u.real[0][3] * amp_pt_i[33] + u.imag[0][3] * amp_pt_r[33]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[32]   + u.imag[1][2] * amp_pt_r[32]) +
                           (u.real[1][3] * amp_pt_i[33] + u.imag[1][3] * amp_pt_r[33]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[32]   + u.imag[2][2] * amp_pt_r[32]) +
                           (u.real[2][3] * amp_pt_i[33] + u.imag[2][3] * amp_pt_r[33]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[32]   + u.imag[3][2] * amp_pt_r[32]) +
                           (u.real[3][3] * amp_pt_i[33] + u.imag[3][3] * amp_pt_r[33]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[32] = amp_r[2];
                amp_pt_r[33] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[32] = amp_i[2];
                amp_pt_i[33] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[65] - u.imag[0][3] * amp_pt_i[65]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[65] - u.imag[1][3] * amp_pt_i[65]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[65] - u.imag[2][3] * amp_pt_i[65]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[65] - u.imag[3][3] * amp_pt_i[65]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[65] + u.imag[0][3] * amp_pt_r[65]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[65] + u.imag[1][3] * amp_pt_r[65]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[65] + u.imag[2][3] * amp_pt_r[65]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[65] + u.imag[3][3] * amp_pt_r[65]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[65] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[65] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[129] - u.imag[0][3] * amp_pt_i[129]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[129] - u.imag[1][3] * amp_pt_i[129]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[129] - u.imag[2][3] * amp_pt_i[129]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[129] - u.imag[3][3] * amp_pt_i[129]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[129] + u.imag[0][3] * amp_pt_r[129]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[129] + u.imag[1][3] * amp_pt_r[129]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[129] + u.imag[2][3] * amp_pt_r[129]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[129] + u.imag[3][3] * amp_pt_r[129]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[129] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[129] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[257] - u.imag[0][3] * amp_pt_i[257]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[257] - u.imag[1][3] * amp_pt_i[257]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[257] - u.imag[2][3] * amp_pt_i[257]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[257] - u.imag[3][3] * amp_pt_i[257]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[257] + u.imag[0][3] * amp_pt_r[257]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[257] + u.imag[1][3] * amp_pt_r[257]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[257] + u.imag[2][3] * amp_pt_r[257]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[257] + u.imag[3][3] * amp_pt_r[257]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[257] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[257] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_0(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (1 << 1))
        {
            for (size_t ttt = 0; ttt < 1; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[1]   - u.imag[0][1] * amp_pt_i[1]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[513] - u.imag[0][3] * amp_pt_i[513]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[1]   - u.imag[1][1] * amp_pt_i[1]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[513] - u.imag[1][3] * amp_pt_i[513]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[1]   - u.imag[2][1] * amp_pt_i[1]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[513] - u.imag[2][3] * amp_pt_i[513]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[1]   - u.imag[3][1] * amp_pt_i[1]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[513] - u.imag[3][3] * amp_pt_i[513]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[1]   + u.imag[0][1] * amp_pt_r[1]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[513] + u.imag[0][3] * amp_pt_r[513]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[1]   + u.imag[1][1] * amp_pt_r[1]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[513] + u.imag[1][3] * amp_pt_r[513]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[1]   + u.imag[2][1] * amp_pt_r[1]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[513] + u.imag[2][3] * amp_pt_r[513]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[1]   + u.imag[3][1] * amp_pt_r[1]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[513] + u.imag[3][3] * amp_pt_r[513]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[1] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[513] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[1] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[513] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_2_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (4 << 1))
    {
        for (size_t tt = 0; tt < 4; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[4]   - u.imag[0][2] * amp_pt_i[4]) +
                           (u.real[0][3] * amp_pt_r[6] - u.imag[0][3] * amp_pt_i[6]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[4]   - u.imag[1][2] * amp_pt_i[4]) +
                           (u.real[1][3] * amp_pt_r[6] - u.imag[1][3] * amp_pt_i[6]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[4]   - u.imag[2][2] * amp_pt_i[4]) +
                           (u.real[2][3] * amp_pt_r[6] - u.imag[2][3] * amp_pt_i[6]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[4]   - u.imag[3][2] * amp_pt_i[4]) +
                           (u.real[3][3] * amp_pt_r[6] - u.imag[3][3] * amp_pt_i[6]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[4]   + u.imag[0][2] * amp_pt_r[4]) +
                           (u.real[0][3] * amp_pt_i[6] + u.imag[0][3] * amp_pt_r[6]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[4]   + u.imag[1][2] * amp_pt_r[4]) +
                           (u.real[1][3] * amp_pt_i[6] + u.imag[1][3] * amp_pt_r[6]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[4]   + u.imag[2][2] * amp_pt_r[4]) +
                           (u.real[2][3] * amp_pt_i[6] + u.imag[2][3] * amp_pt_r[6]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[4]   + u.imag[3][2] * amp_pt_r[4]) +
                           (u.real[3][3] * amp_pt_i[6] + u.imag[3][3] * amp_pt_r[6]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[4] = amp_r[2];
                amp_pt_r[6] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[4] = amp_i[2];
                amp_pt_i[6] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_3_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (8 << 1))
    {
        for (size_t tt = 0; tt < 8; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[8]   - u.imag[0][2] * amp_pt_i[8]) +
                           (u.real[0][3] * amp_pt_r[10] - u.imag[0][3] * amp_pt_i[10]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[8]   - u.imag[1][2] * amp_pt_i[8]) +
                           (u.real[1][3] * amp_pt_r[10] - u.imag[1][3] * amp_pt_i[10]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[8]   - u.imag[2][2] * amp_pt_i[8]) +
                           (u.real[2][3] * amp_pt_r[10] - u.imag[2][3] * amp_pt_i[10]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[8]   - u.imag[3][2] * amp_pt_i[8]) +
                           (u.real[3][3] * amp_pt_r[10] - u.imag[3][3] * amp_pt_i[10]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[8]   + u.imag[0][2] * amp_pt_r[8]) +
                           (u.real[0][3] * amp_pt_i[10] + u.imag[0][3] * amp_pt_r[10]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[8]   + u.imag[1][2] * amp_pt_r[8]) +
                           (u.real[1][3] * amp_pt_i[10] + u.imag[1][3] * amp_pt_r[10]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[8]   + u.imag[2][2] * amp_pt_r[8]) +
                           (u.real[2][3] * amp_pt_i[10] + u.imag[2][3] * amp_pt_r[10]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[8]   + u.imag[3][2] * amp_pt_r[8]) +
                           (u.real[3][3] * amp_pt_i[10] + u.imag[3][3] * amp_pt_r[10]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[8] = amp_r[2];
                amp_pt_r[10] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[8] = amp_i[2];
                amp_pt_i[10] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_4_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (16 << 1))
    {
        for (size_t tt = 0; tt < 16; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[16]   - u.imag[0][2] * amp_pt_i[16]) +
                           (u.real[0][3] * amp_pt_r[18] - u.imag[0][3] * amp_pt_i[18]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[16]   - u.imag[1][2] * amp_pt_i[16]) +
                           (u.real[1][3] * amp_pt_r[18] - u.imag[1][3] * amp_pt_i[18]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[16]   - u.imag[2][2] * amp_pt_i[16]) +
                           (u.real[2][3] * amp_pt_r[18] - u.imag[2][3] * amp_pt_i[18]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[16]   - u.imag[3][2] * amp_pt_i[16]) +
                           (u.real[3][3] * amp_pt_r[18] - u.imag[3][3] * amp_pt_i[18]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[16]   + u.imag[0][2] * amp_pt_r[16]) +
                           (u.real[0][3] * amp_pt_i[18] + u.imag[0][3] * amp_pt_r[18]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[16]   + u.imag[1][2] * amp_pt_r[16]) +
                           (u.real[1][3] * amp_pt_i[18] + u.imag[1][3] * amp_pt_r[18]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[16]   + u.imag[2][2] * amp_pt_r[16]) +
                           (u.real[2][3] * amp_pt_i[18] + u.imag[2][3] * amp_pt_r[18]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[16]   + u.imag[3][2] * amp_pt_r[16]) +
                           (u.real[3][3] * amp_pt_i[18] + u.imag[3][3] * amp_pt_r[18]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[16] = amp_r[2];
                amp_pt_r[18] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[16] = amp_i[2];
                amp_pt_i[18] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_5_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (32 << 1))
    {
        for (size_t tt = 0; tt < 32; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[32]   - u.imag[0][2] * amp_pt_i[32]) +
                           (u.real[0][3] * amp_pt_r[34] - u.imag[0][3] * amp_pt_i[34]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[32]   - u.imag[1][2] * amp_pt_i[32]) +
                           (u.real[1][3] * amp_pt_r[34] - u.imag[1][3] * amp_pt_i[34]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[32]   - u.imag[2][2] * amp_pt_i[32]) +
                           (u.real[2][3] * amp_pt_r[34] - u.imag[2][3] * amp_pt_i[34]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[32]   - u.imag[3][2] * amp_pt_i[32]) +
                           (u.real[3][3] * amp_pt_r[34] - u.imag[3][3] * amp_pt_i[34]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[32]   + u.imag[0][2] * amp_pt_r[32]) +
                           (u.real[0][3] * amp_pt_i[34] + u.imag[0][3] * amp_pt_r[34]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[32]   + u.imag[1][2] * amp_pt_r[32]) +
                           (u.real[1][3] * amp_pt_i[34] + u.imag[1][3] * amp_pt_r[34]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[32]   + u.imag[2][2] * amp_pt_r[32]) +
                           (u.real[2][3] * amp_pt_i[34] + u.imag[2][3] * amp_pt_r[34]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[32]   + u.imag[3][2] * amp_pt_r[32]) +
                           (u.real[3][3] * amp_pt_i[34] + u.imag[3][3] * amp_pt_r[34]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[32] = amp_r[2];
                amp_pt_r[34] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[32] = amp_i[2];
                amp_pt_i[34] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[66] - u.imag[0][3] * amp_pt_i[66]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[66] - u.imag[1][3] * amp_pt_i[66]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[66] - u.imag[2][3] * amp_pt_i[66]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[66] - u.imag[3][3] * amp_pt_i[66]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[66] + u.imag[0][3] * amp_pt_r[66]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[66] + u.imag[1][3] * amp_pt_r[66]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[66] + u.imag[2][3] * amp_pt_r[66]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[66] + u.imag[3][3] * amp_pt_r[66]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[66] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[66] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[130] - u.imag[0][3] * amp_pt_i[130]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[130] - u.imag[1][3] * amp_pt_i[130]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[130] - u.imag[2][3] * amp_pt_i[130]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[130] - u.imag[3][3] * amp_pt_i[130]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[130] + u.imag[0][3] * amp_pt_r[130]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[130] + u.imag[1][3] * amp_pt_r[130]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[130] + u.imag[2][3] * amp_pt_r[130]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[130] + u.imag[3][3] * amp_pt_r[130]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[130] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[130] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[258] - u.imag[0][3] * amp_pt_i[258]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[258] - u.imag[1][3] * amp_pt_i[258]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[258] - u.imag[2][3] * amp_pt_i[258]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[258] - u.imag[3][3] * amp_pt_i[258]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[258] + u.imag[0][3] * amp_pt_r[258]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[258] + u.imag[1][3] * amp_pt_r[258]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[258] + u.imag[2][3] * amp_pt_r[258]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[258] + u.imag[3][3] * amp_pt_r[258]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[258] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[258] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_1(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (2 << 1))
        {
            for (size_t ttt = 0; ttt < 2; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[2]   - u.imag[0][1] * amp_pt_i[2]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[514] - u.imag[0][3] * amp_pt_i[514]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[2]   - u.imag[1][1] * amp_pt_i[2]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[514] - u.imag[1][3] * amp_pt_i[514]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[2]   - u.imag[2][1] * amp_pt_i[2]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[514] - u.imag[2][3] * amp_pt_i[514]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[2]   - u.imag[3][1] * amp_pt_i[2]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[514] - u.imag[3][3] * amp_pt_i[514]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[2]   + u.imag[0][1] * amp_pt_r[2]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[514] + u.imag[0][3] * amp_pt_r[514]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[2]   + u.imag[1][1] * amp_pt_r[2]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[514] + u.imag[1][3] * amp_pt_r[514]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[2]   + u.imag[2][1] * amp_pt_r[2]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[514] + u.imag[2][3] * amp_pt_r[514]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[2]   + u.imag[3][1] * amp_pt_r[2]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[514] + u.imag[3][3] * amp_pt_r[514]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[2] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[514] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[2] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[514] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_3_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (8 << 1))
    {
        for (size_t tt = 0; tt < 8; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[8]   - u.imag[0][2] * amp_pt_i[8]) +
                           (u.real[0][3] * amp_pt_r[12] - u.imag[0][3] * amp_pt_i[12]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[8]   - u.imag[1][2] * amp_pt_i[8]) +
                           (u.real[1][3] * amp_pt_r[12] - u.imag[1][3] * amp_pt_i[12]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[8]   - u.imag[2][2] * amp_pt_i[8]) +
                           (u.real[2][3] * amp_pt_r[12] - u.imag[2][3] * amp_pt_i[12]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[8]   - u.imag[3][2] * amp_pt_i[8]) +
                           (u.real[3][3] * amp_pt_r[12] - u.imag[3][3] * amp_pt_i[12]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[8]   + u.imag[0][2] * amp_pt_r[8]) +
                           (u.real[0][3] * amp_pt_i[12] + u.imag[0][3] * amp_pt_r[12]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[8]   + u.imag[1][2] * amp_pt_r[8]) +
                           (u.real[1][3] * amp_pt_i[12] + u.imag[1][3] * amp_pt_r[12]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[8]   + u.imag[2][2] * amp_pt_r[8]) +
                           (u.real[2][3] * amp_pt_i[12] + u.imag[2][3] * amp_pt_r[12]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[8]   + u.imag[3][2] * amp_pt_r[8]) +
                           (u.real[3][3] * amp_pt_i[12] + u.imag[3][3] * amp_pt_r[12]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[8] = amp_r[2];
                amp_pt_r[12] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[8] = amp_i[2];
                amp_pt_i[12] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_4_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (16 << 1))
    {
        for (size_t tt = 0; tt < 16; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[16]   - u.imag[0][2] * amp_pt_i[16]) +
                           (u.real[0][3] * amp_pt_r[20] - u.imag[0][3] * amp_pt_i[20]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[16]   - u.imag[1][2] * amp_pt_i[16]) +
                           (u.real[1][3] * amp_pt_r[20] - u.imag[1][3] * amp_pt_i[20]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[16]   - u.imag[2][2] * amp_pt_i[16]) +
                           (u.real[2][3] * amp_pt_r[20] - u.imag[2][3] * amp_pt_i[20]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[16]   - u.imag[3][2] * amp_pt_i[16]) +
                           (u.real[3][3] * amp_pt_r[20] - u.imag[3][3] * amp_pt_i[20]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[16]   + u.imag[0][2] * amp_pt_r[16]) +
                           (u.real[0][3] * amp_pt_i[20] + u.imag[0][3] * amp_pt_r[20]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[16]   + u.imag[1][2] * amp_pt_r[16]) +
                           (u.real[1][3] * amp_pt_i[20] + u.imag[1][3] * amp_pt_r[20]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[16]   + u.imag[2][2] * amp_pt_r[16]) +
                           (u.real[2][3] * amp_pt_i[20] + u.imag[2][3] * amp_pt_r[20]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[16]   + u.imag[3][2] * amp_pt_r[16]) +
                           (u.real[3][3] * amp_pt_i[20] + u.imag[3][3] * amp_pt_r[20]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[16] = amp_r[2];
                amp_pt_r[20] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[16] = amp_i[2];
                amp_pt_i[20] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_5_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (32 << 1))
    {
        for (size_t tt = 0; tt < 32; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[32]   - u.imag[0][2] * amp_pt_i[32]) +
                           (u.real[0][3] * amp_pt_r[36] - u.imag[0][3] * amp_pt_i[36]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[32]   - u.imag[1][2] * amp_pt_i[32]) +
                           (u.real[1][3] * amp_pt_r[36] - u.imag[1][3] * amp_pt_i[36]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[32]   - u.imag[2][2] * amp_pt_i[32]) +
                           (u.real[2][3] * amp_pt_r[36] - u.imag[2][3] * amp_pt_i[36]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[32]   - u.imag[3][2] * amp_pt_i[32]) +
                           (u.real[3][3] * amp_pt_r[36] - u.imag[3][3] * amp_pt_i[36]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[32]   + u.imag[0][2] * amp_pt_r[32]) +
                           (u.real[0][3] * amp_pt_i[36] + u.imag[0][3] * amp_pt_r[36]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[32]   + u.imag[1][2] * amp_pt_r[32]) +
                           (u.real[1][3] * amp_pt_i[36] + u.imag[1][3] * amp_pt_r[36]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[32]   + u.imag[2][2] * amp_pt_r[32]) +
                           (u.real[2][3] * amp_pt_i[36] + u.imag[2][3] * amp_pt_r[36]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[32]   + u.imag[3][2] * amp_pt_r[32]) +
                           (u.real[3][3] * amp_pt_i[36] + u.imag[3][3] * amp_pt_r[36]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[32] = amp_r[2];
                amp_pt_r[36] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[32] = amp_i[2];
                amp_pt_i[36] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[68] - u.imag[0][3] * amp_pt_i[68]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[68] - u.imag[1][3] * amp_pt_i[68]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[68] - u.imag[2][3] * amp_pt_i[68]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[68] - u.imag[3][3] * amp_pt_i[68]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[68] + u.imag[0][3] * amp_pt_r[68]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[68] + u.imag[1][3] * amp_pt_r[68]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[68] + u.imag[2][3] * amp_pt_r[68]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[68] + u.imag[3][3] * amp_pt_r[68]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[68] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[68] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[132] - u.imag[0][3] * amp_pt_i[132]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[132] - u.imag[1][3] * amp_pt_i[132]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[132] - u.imag[2][3] * amp_pt_i[132]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[132] - u.imag[3][3] * amp_pt_i[132]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[132] + u.imag[0][3] * amp_pt_r[132]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[132] + u.imag[1][3] * amp_pt_r[132]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[132] + u.imag[2][3] * amp_pt_r[132]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[132] + u.imag[3][3] * amp_pt_r[132]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[132] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[132] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[260] - u.imag[0][3] * amp_pt_i[260]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[260] - u.imag[1][3] * amp_pt_i[260]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[260] - u.imag[2][3] * amp_pt_i[260]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[260] - u.imag[3][3] * amp_pt_i[260]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[260] + u.imag[0][3] * amp_pt_r[260]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[260] + u.imag[1][3] * amp_pt_r[260]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[260] + u.imag[2][3] * amp_pt_r[260]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[260] + u.imag[3][3] * amp_pt_r[260]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[260] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[260] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_2(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (4 << 1))
        {
            for (size_t ttt = 0; ttt < 4; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[4]   - u.imag[0][1] * amp_pt_i[4]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[516] - u.imag[0][3] * amp_pt_i[516]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[4]   - u.imag[1][1] * amp_pt_i[4]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[516] - u.imag[1][3] * amp_pt_i[516]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[4]   - u.imag[2][1] * amp_pt_i[4]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[516] - u.imag[2][3] * amp_pt_i[516]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[4]   - u.imag[3][1] * amp_pt_i[4]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[516] - u.imag[3][3] * amp_pt_i[516]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[4]   + u.imag[0][1] * amp_pt_r[4]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[516] + u.imag[0][3] * amp_pt_r[516]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[4]   + u.imag[1][1] * amp_pt_r[4]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[516] + u.imag[1][3] * amp_pt_r[516]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[4]   + u.imag[2][1] * amp_pt_r[4]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[516] + u.imag[2][3] * amp_pt_r[516]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[4]   + u.imag[3][1] * amp_pt_r[4]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[516] + u.imag[3][3] * amp_pt_r[516]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[4] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[516] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[4] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[516] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_4_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (16 << 1))
    {
        for (size_t tt = 0; tt < 16; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[16]   - u.imag[0][2] * amp_pt_i[16]) +
                           (u.real[0][3] * amp_pt_r[24] - u.imag[0][3] * amp_pt_i[24]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[16]   - u.imag[1][2] * amp_pt_i[16]) +
                           (u.real[1][3] * amp_pt_r[24] - u.imag[1][3] * amp_pt_i[24]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[16]   - u.imag[2][2] * amp_pt_i[16]) +
                           (u.real[2][3] * amp_pt_r[24] - u.imag[2][3] * amp_pt_i[24]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[16]   - u.imag[3][2] * amp_pt_i[16]) +
                           (u.real[3][3] * amp_pt_r[24] - u.imag[3][3] * amp_pt_i[24]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[16]   + u.imag[0][2] * amp_pt_r[16]) +
                           (u.real[0][3] * amp_pt_i[24] + u.imag[0][3] * amp_pt_r[24]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[16]   + u.imag[1][2] * amp_pt_r[16]) +
                           (u.real[1][3] * amp_pt_i[24] + u.imag[1][3] * amp_pt_r[24]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[16]   + u.imag[2][2] * amp_pt_r[16]) +
                           (u.real[2][3] * amp_pt_i[24] + u.imag[2][3] * amp_pt_r[24]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[16]   + u.imag[3][2] * amp_pt_r[16]) +
                           (u.real[3][3] * amp_pt_i[24] + u.imag[3][3] * amp_pt_r[24]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[16] = amp_r[2];
                amp_pt_r[24] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[16] = amp_i[2];
                amp_pt_i[24] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_5_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (32 << 1))
    {
        for (size_t tt = 0; tt < 32; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[32]   - u.imag[0][2] * amp_pt_i[32]) +
                           (u.real[0][3] * amp_pt_r[40] - u.imag[0][3] * amp_pt_i[40]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[32]   - u.imag[1][2] * amp_pt_i[32]) +
                           (u.real[1][3] * amp_pt_r[40] - u.imag[1][3] * amp_pt_i[40]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[32]   - u.imag[2][2] * amp_pt_i[32]) +
                           (u.real[2][3] * amp_pt_r[40] - u.imag[2][3] * amp_pt_i[40]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[32]   - u.imag[3][2] * amp_pt_i[32]) +
                           (u.real[3][3] * amp_pt_r[40] - u.imag[3][3] * amp_pt_i[40]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[32]   + u.imag[0][2] * amp_pt_r[32]) +
                           (u.real[0][3] * amp_pt_i[40] + u.imag[0][3] * amp_pt_r[40]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[32]   + u.imag[1][2] * amp_pt_r[32]) +
                           (u.real[1][3] * amp_pt_i[40] + u.imag[1][3] * amp_pt_r[40]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[32]   + u.imag[2][2] * amp_pt_r[32]) +
                           (u.real[2][3] * amp_pt_i[40] + u.imag[2][3] * amp_pt_r[40]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[32]   + u.imag[3][2] * amp_pt_r[32]) +
                           (u.real[3][3] * amp_pt_i[40] + u.imag[3][3] * amp_pt_r[40]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[32] = amp_r[2];
                amp_pt_r[40] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[32] = amp_i[2];
                amp_pt_i[40] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[72] - u.imag[0][3] * amp_pt_i[72]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[72] - u.imag[1][3] * amp_pt_i[72]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[72] - u.imag[2][3] * amp_pt_i[72]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[72] - u.imag[3][3] * amp_pt_i[72]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[72] + u.imag[0][3] * amp_pt_r[72]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[72] + u.imag[1][3] * amp_pt_r[72]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[72] + u.imag[2][3] * amp_pt_r[72]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[72] + u.imag[3][3] * amp_pt_r[72]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[72] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[72] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[136] - u.imag[0][3] * amp_pt_i[136]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[136] - u.imag[1][3] * amp_pt_i[136]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[136] - u.imag[2][3] * amp_pt_i[136]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[136] - u.imag[3][3] * amp_pt_i[136]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[136] + u.imag[0][3] * amp_pt_r[136]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[136] + u.imag[1][3] * amp_pt_r[136]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[136] + u.imag[2][3] * amp_pt_r[136]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[136] + u.imag[3][3] * amp_pt_r[136]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[136] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[136] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[264] - u.imag[0][3] * amp_pt_i[264]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[264] - u.imag[1][3] * amp_pt_i[264]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[264] - u.imag[2][3] * amp_pt_i[264]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[264] - u.imag[3][3] * amp_pt_i[264]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[264] + u.imag[0][3] * amp_pt_r[264]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[264] + u.imag[1][3] * amp_pt_r[264]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[264] + u.imag[2][3] * amp_pt_r[264]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[264] + u.imag[3][3] * amp_pt_r[264]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[264] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[264] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_3(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (8 << 1))
        {
            for (size_t ttt = 0; ttt < 8; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[8]   - u.imag[0][1] * amp_pt_i[8]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[520] - u.imag[0][3] * amp_pt_i[520]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[8]   - u.imag[1][1] * amp_pt_i[8]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[520] - u.imag[1][3] * amp_pt_i[520]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[8]   - u.imag[2][1] * amp_pt_i[8]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[520] - u.imag[2][3] * amp_pt_i[520]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[8]   - u.imag[3][1] * amp_pt_i[8]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[520] - u.imag[3][3] * amp_pt_i[520]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[8]   + u.imag[0][1] * amp_pt_r[8]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[520] + u.imag[0][3] * amp_pt_r[520]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[8]   + u.imag[1][1] * amp_pt_r[8]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[520] + u.imag[1][3] * amp_pt_r[520]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[8]   + u.imag[2][1] * amp_pt_r[8]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[520] + u.imag[2][3] * amp_pt_r[520]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[8]   + u.imag[3][1] * amp_pt_r[8]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[520] + u.imag[3][3] * amp_pt_r[520]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[8] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[520] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[8] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[520] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_5_4(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (32 << 1))
    {
        for (size_t tt = 0; tt < 32; tt += (16 << 1))
        {
            for (size_t ttt = 0; ttt < 16; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[16]   - u.imag[0][1] * amp_pt_i[16]) +
                           (u.real[0][2] * amp_pt_r[32]   - u.imag[0][2] * amp_pt_i[32]) +
                           (u.real[0][3] * amp_pt_r[48] - u.imag[0][3] * amp_pt_i[48]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[16]   - u.imag[1][1] * amp_pt_i[16]) +
                           (u.real[1][2] * amp_pt_r[32]   - u.imag[1][2] * amp_pt_i[32]) +
                           (u.real[1][3] * amp_pt_r[48] - u.imag[1][3] * amp_pt_i[48]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[16]   - u.imag[2][1] * amp_pt_i[16]) +
                           (u.real[2][2] * amp_pt_r[32]   - u.imag[2][2] * amp_pt_i[32]) +
                           (u.real[2][3] * amp_pt_r[48] - u.imag[2][3] * amp_pt_i[48]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[16]   - u.imag[3][1] * amp_pt_i[16]) +
                           (u.real[3][2] * amp_pt_r[32]   - u.imag[3][2] * amp_pt_i[32]) +
                           (u.real[3][3] * amp_pt_r[48] - u.imag[3][3] * amp_pt_i[48]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[16]   + u.imag[0][1] * amp_pt_r[16]) +
                           (u.real[0][2] * amp_pt_i[32]   + u.imag[0][2] * amp_pt_r[32]) +
                           (u.real[0][3] * amp_pt_i[48] + u.imag[0][3] * amp_pt_r[48]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[16]   + u.imag[1][1] * amp_pt_r[16]) +
                           (u.real[1][2] * amp_pt_i[32]   + u.imag[1][2] * amp_pt_r[32]) +
                           (u.real[1][3] * amp_pt_i[48] + u.imag[1][3] * amp_pt_r[48]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[16]   + u.imag[2][1] * amp_pt_r[16]) +
                           (u.real[2][2] * amp_pt_i[32]   + u.imag[2][2] * amp_pt_r[32]) +
                           (u.real[2][3] * amp_pt_i[48] + u.imag[2][3] * amp_pt_r[48]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[16]   + u.imag[3][1] * amp_pt_r[16]) +
                           (u.real[3][2] * amp_pt_i[32]   + u.imag[3][2] * amp_pt_r[32]) +
                           (u.real[3][3] * amp_pt_i[48] + u.imag[3][3] * amp_pt_r[48]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[16] = amp_r[1];
                amp_pt_r[32] = amp_r[2];
                amp_pt_r[48] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[16] = amp_i[1];
                amp_pt_i[32] = amp_i[2];
                amp_pt_i[48] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_4(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (16 << 1))
        {
            for (size_t ttt = 0; ttt < 16; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[16]   - u.imag[0][1] * amp_pt_i[16]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[80] - u.imag[0][3] * amp_pt_i[80]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[16]   - u.imag[1][1] * amp_pt_i[16]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[80] - u.imag[1][3] * amp_pt_i[80]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[16]   - u.imag[2][1] * amp_pt_i[16]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[80] - u.imag[2][3] * amp_pt_i[80]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[16]   - u.imag[3][1] * amp_pt_i[16]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[80] - u.imag[3][3] * amp_pt_i[80]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[16]   + u.imag[0][1] * amp_pt_r[16]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[80] + u.imag[0][3] * amp_pt_r[80]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[16]   + u.imag[1][1] * amp_pt_r[16]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[80] + u.imag[1][3] * amp_pt_r[80]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[16]   + u.imag[2][1] * amp_pt_r[16]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[80] + u.imag[2][3] * amp_pt_r[80]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[16]   + u.imag[3][1] * amp_pt_r[16]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[80] + u.imag[3][3] * amp_pt_r[80]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[16] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[80] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[16] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[80] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_4(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (16 << 1))
        {
            for (size_t ttt = 0; ttt < 16; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[16]   - u.imag[0][1] * amp_pt_i[16]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[144] - u.imag[0][3] * amp_pt_i[144]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[16]   - u.imag[1][1] * amp_pt_i[16]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[144] - u.imag[1][3] * amp_pt_i[144]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[16]   - u.imag[2][1] * amp_pt_i[16]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[144] - u.imag[2][3] * amp_pt_i[144]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[16]   - u.imag[3][1] * amp_pt_i[16]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[144] - u.imag[3][3] * amp_pt_i[144]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[16]   + u.imag[0][1] * amp_pt_r[16]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[144] + u.imag[0][3] * amp_pt_r[144]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[16]   + u.imag[1][1] * amp_pt_r[16]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[144] + u.imag[1][3] * amp_pt_r[144]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[16]   + u.imag[2][1] * amp_pt_r[16]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[144] + u.imag[2][3] * amp_pt_r[144]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[16]   + u.imag[3][1] * amp_pt_r[16]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[144] + u.imag[3][3] * amp_pt_r[144]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[16] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[144] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[16] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[144] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_4(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (16 << 1))
        {
            for (size_t ttt = 0; ttt < 16; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[16]   - u.imag[0][1] * amp_pt_i[16]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[272] - u.imag[0][3] * amp_pt_i[272]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[16]   - u.imag[1][1] * amp_pt_i[16]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[272] - u.imag[1][3] * amp_pt_i[272]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[16]   - u.imag[2][1] * amp_pt_i[16]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[272] - u.imag[2][3] * amp_pt_i[272]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[16]   - u.imag[3][1] * amp_pt_i[16]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[272] - u.imag[3][3] * amp_pt_i[272]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[16]   + u.imag[0][1] * amp_pt_r[16]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[272] + u.imag[0][3] * amp_pt_r[272]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[16]   + u.imag[1][1] * amp_pt_r[16]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[272] + u.imag[1][3] * amp_pt_r[272]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[16]   + u.imag[2][1] * amp_pt_r[16]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[272] + u.imag[2][3] * amp_pt_r[272]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[16]   + u.imag[3][1] * amp_pt_r[16]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[272] + u.imag[3][3] * amp_pt_r[272]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[16] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[272] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[16] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[272] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_4(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (16 << 1))
        {
            for (size_t ttt = 0; ttt < 16; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[16]   - u.imag[0][1] * amp_pt_i[16]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[528] - u.imag[0][3] * amp_pt_i[528]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[16]   - u.imag[1][1] * amp_pt_i[16]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[528] - u.imag[1][3] * amp_pt_i[528]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[16]   - u.imag[2][1] * amp_pt_i[16]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[528] - u.imag[2][3] * amp_pt_i[528]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[16]   - u.imag[3][1] * amp_pt_i[16]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[528] - u.imag[3][3] * amp_pt_i[528]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[16]   + u.imag[0][1] * amp_pt_r[16]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[528] + u.imag[0][3] * amp_pt_r[528]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[16]   + u.imag[1][1] * amp_pt_r[16]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[528] + u.imag[1][3] * amp_pt_r[528]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[16]   + u.imag[2][1] * amp_pt_r[16]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[528] + u.imag[2][3] * amp_pt_r[528]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[16]   + u.imag[3][1] * amp_pt_r[16]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[528] + u.imag[3][3] * amp_pt_r[528]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[16] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[528] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[16] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[528] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_6_5(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (64 << 1))
    {
        for (size_t tt = 0; tt < 64; tt += (32 << 1))
        {
            for (size_t ttt = 0; ttt < 32; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[32]   - u.imag[0][1] * amp_pt_i[32]) +
                           (u.real[0][2] * amp_pt_r[64]   - u.imag[0][2] * amp_pt_i[64]) +
                           (u.real[0][3] * amp_pt_r[96] - u.imag[0][3] * amp_pt_i[96]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[32]   - u.imag[1][1] * amp_pt_i[32]) +
                           (u.real[1][2] * amp_pt_r[64]   - u.imag[1][2] * amp_pt_i[64]) +
                           (u.real[1][3] * amp_pt_r[96] - u.imag[1][3] * amp_pt_i[96]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[32]   - u.imag[2][1] * amp_pt_i[32]) +
                           (u.real[2][2] * amp_pt_r[64]   - u.imag[2][2] * amp_pt_i[64]) +
                           (u.real[2][3] * amp_pt_r[96] - u.imag[2][3] * amp_pt_i[96]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[32]   - u.imag[3][1] * amp_pt_i[32]) +
                           (u.real[3][2] * amp_pt_r[64]   - u.imag[3][2] * amp_pt_i[64]) +
                           (u.real[3][3] * amp_pt_r[96] - u.imag[3][3] * amp_pt_i[96]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[32]   + u.imag[0][1] * amp_pt_r[32]) +
                           (u.real[0][2] * amp_pt_i[64]   + u.imag[0][2] * amp_pt_r[64]) +
                           (u.real[0][3] * amp_pt_i[96] + u.imag[0][3] * amp_pt_r[96]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[32]   + u.imag[1][1] * amp_pt_r[32]) +
                           (u.real[1][2] * amp_pt_i[64]   + u.imag[1][2] * amp_pt_r[64]) +
                           (u.real[1][3] * amp_pt_i[96] + u.imag[1][3] * amp_pt_r[96]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[32]   + u.imag[2][1] * amp_pt_r[32]) +
                           (u.real[2][2] * amp_pt_i[64]   + u.imag[2][2] * amp_pt_r[64]) +
                           (u.real[2][3] * amp_pt_i[96] + u.imag[2][3] * amp_pt_r[96]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[32]   + u.imag[3][1] * amp_pt_r[32]) +
                           (u.real[3][2] * amp_pt_i[64]   + u.imag[3][2] * amp_pt_r[64]) +
                           (u.real[3][3] * amp_pt_i[96] + u.imag[3][3] * amp_pt_r[96]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[32] = amp_r[1];
                amp_pt_r[64] = amp_r[2];
                amp_pt_r[96] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[32] = amp_i[1];
                amp_pt_i[64] = amp_i[2];
                amp_pt_i[96] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_5(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (32 << 1))
        {
            for (size_t ttt = 0; ttt < 32; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[32]   - u.imag[0][1] * amp_pt_i[32]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[160] - u.imag[0][3] * amp_pt_i[160]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[32]   - u.imag[1][1] * amp_pt_i[32]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[160] - u.imag[1][3] * amp_pt_i[160]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[32]   - u.imag[2][1] * amp_pt_i[32]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[160] - u.imag[2][3] * amp_pt_i[160]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[32]   - u.imag[3][1] * amp_pt_i[32]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[160] - u.imag[3][3] * amp_pt_i[160]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[32]   + u.imag[0][1] * amp_pt_r[32]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[160] + u.imag[0][3] * amp_pt_r[160]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[32]   + u.imag[1][1] * amp_pt_r[32]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[160] + u.imag[1][3] * amp_pt_r[160]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[32]   + u.imag[2][1] * amp_pt_r[32]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[160] + u.imag[2][3] * amp_pt_r[160]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[32]   + u.imag[3][1] * amp_pt_r[32]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[160] + u.imag[3][3] * amp_pt_r[160]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[32] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[160] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[32] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[160] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_5(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (32 << 1))
        {
            for (size_t ttt = 0; ttt < 32; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[32]   - u.imag[0][1] * amp_pt_i[32]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[288] - u.imag[0][3] * amp_pt_i[288]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[32]   - u.imag[1][1] * amp_pt_i[32]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[288] - u.imag[1][3] * amp_pt_i[288]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[32]   - u.imag[2][1] * amp_pt_i[32]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[288] - u.imag[2][3] * amp_pt_i[288]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[32]   - u.imag[3][1] * amp_pt_i[32]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[288] - u.imag[3][3] * amp_pt_i[288]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[32]   + u.imag[0][1] * amp_pt_r[32]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[288] + u.imag[0][3] * amp_pt_r[288]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[32]   + u.imag[1][1] * amp_pt_r[32]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[288] + u.imag[1][3] * amp_pt_r[288]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[32]   + u.imag[2][1] * amp_pt_r[32]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[288] + u.imag[2][3] * amp_pt_r[288]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[32]   + u.imag[3][1] * amp_pt_r[32]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[288] + u.imag[3][3] * amp_pt_r[288]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[32] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[288] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[32] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[288] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_5(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (32 << 1))
        {
            for (size_t ttt = 0; ttt < 32; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[32]   - u.imag[0][1] * amp_pt_i[32]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[544] - u.imag[0][3] * amp_pt_i[544]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[32]   - u.imag[1][1] * amp_pt_i[32]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[544] - u.imag[1][3] * amp_pt_i[544]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[32]   - u.imag[2][1] * amp_pt_i[32]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[544] - u.imag[2][3] * amp_pt_i[544]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[32]   - u.imag[3][1] * amp_pt_i[32]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[544] - u.imag[3][3] * amp_pt_i[544]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[32]   + u.imag[0][1] * amp_pt_r[32]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[544] + u.imag[0][3] * amp_pt_r[544]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[32]   + u.imag[1][1] * amp_pt_r[32]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[544] + u.imag[1][3] * amp_pt_r[544]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[32]   + u.imag[2][1] * amp_pt_r[32]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[544] + u.imag[2][3] * amp_pt_r[544]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[32]   + u.imag[3][1] * amp_pt_r[32]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[544] + u.imag[3][3] * amp_pt_r[544]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[32] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[544] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[32] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[544] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_7_6(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (128 << 1))
    {
        for (size_t tt = 0; tt < 128; tt += (64 << 1))
        {
            for (size_t ttt = 0; ttt < 64; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[64]   - u.imag[0][1] * amp_pt_i[64]) +
                           (u.real[0][2] * amp_pt_r[128]   - u.imag[0][2] * amp_pt_i[128]) +
                           (u.real[0][3] * amp_pt_r[192] - u.imag[0][3] * amp_pt_i[192]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[64]   - u.imag[1][1] * amp_pt_i[64]) +
                           (u.real[1][2] * amp_pt_r[128]   - u.imag[1][2] * amp_pt_i[128]) +
                           (u.real[1][3] * amp_pt_r[192] - u.imag[1][3] * amp_pt_i[192]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[64]   - u.imag[2][1] * amp_pt_i[64]) +
                           (u.real[2][2] * amp_pt_r[128]   - u.imag[2][2] * amp_pt_i[128]) +
                           (u.real[2][3] * amp_pt_r[192] - u.imag[2][3] * amp_pt_i[192]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[64]   - u.imag[3][1] * amp_pt_i[64]) +
                           (u.real[3][2] * amp_pt_r[128]   - u.imag[3][2] * amp_pt_i[128]) +
                           (u.real[3][3] * amp_pt_r[192] - u.imag[3][3] * amp_pt_i[192]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[64]   + u.imag[0][1] * amp_pt_r[64]) +
                           (u.real[0][2] * amp_pt_i[128]   + u.imag[0][2] * amp_pt_r[128]) +
                           (u.real[0][3] * amp_pt_i[192] + u.imag[0][3] * amp_pt_r[192]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[64]   + u.imag[1][1] * amp_pt_r[64]) +
                           (u.real[1][2] * amp_pt_i[128]   + u.imag[1][2] * amp_pt_r[128]) +
                           (u.real[1][3] * amp_pt_i[192] + u.imag[1][3] * amp_pt_r[192]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[64]   + u.imag[2][1] * amp_pt_r[64]) +
                           (u.real[2][2] * amp_pt_i[128]   + u.imag[2][2] * amp_pt_r[128]) +
                           (u.real[2][3] * amp_pt_i[192] + u.imag[2][3] * amp_pt_r[192]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[64]   + u.imag[3][1] * amp_pt_r[64]) +
                           (u.real[3][2] * amp_pt_i[128]   + u.imag[3][2] * amp_pt_r[128]) +
                           (u.real[3][3] * amp_pt_i[192] + u.imag[3][3] * amp_pt_r[192]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[64] = amp_r[1];
                amp_pt_r[128] = amp_r[2];
                amp_pt_r[192] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[64] = amp_i[1];
                amp_pt_i[128] = amp_i[2];
                amp_pt_i[192] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_6(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (64 << 1))
        {
            for (size_t ttt = 0; ttt < 64; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[64]   - u.imag[0][1] * amp_pt_i[64]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[320] - u.imag[0][3] * amp_pt_i[320]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[64]   - u.imag[1][1] * amp_pt_i[64]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[320] - u.imag[1][3] * amp_pt_i[320]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[64]   - u.imag[2][1] * amp_pt_i[64]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[320] - u.imag[2][3] * amp_pt_i[320]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[64]   - u.imag[3][1] * amp_pt_i[64]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[320] - u.imag[3][3] * amp_pt_i[320]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[64]   + u.imag[0][1] * amp_pt_r[64]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[320] + u.imag[0][3] * amp_pt_r[320]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[64]   + u.imag[1][1] * amp_pt_r[64]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[320] + u.imag[1][3] * amp_pt_r[320]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[64]   + u.imag[2][1] * amp_pt_r[64]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[320] + u.imag[2][3] * amp_pt_r[320]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[64]   + u.imag[3][1] * amp_pt_r[64]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[320] + u.imag[3][3] * amp_pt_r[320]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[64] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[320] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[64] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[320] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_6(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (64 << 1))
        {
            for (size_t ttt = 0; ttt < 64; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[64]   - u.imag[0][1] * amp_pt_i[64]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[576] - u.imag[0][3] * amp_pt_i[576]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[64]   - u.imag[1][1] * amp_pt_i[64]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[576] - u.imag[1][3] * amp_pt_i[576]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[64]   - u.imag[2][1] * amp_pt_i[64]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[576] - u.imag[2][3] * amp_pt_i[576]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[64]   - u.imag[3][1] * amp_pt_i[64]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[576] - u.imag[3][3] * amp_pt_i[576]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[64]   + u.imag[0][1] * amp_pt_r[64]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[576] + u.imag[0][3] * amp_pt_r[576]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[64]   + u.imag[1][1] * amp_pt_r[64]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[576] + u.imag[1][3] * amp_pt_r[576]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[64]   + u.imag[2][1] * amp_pt_r[64]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[576] + u.imag[2][3] * amp_pt_r[576]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[64]   + u.imag[3][1] * amp_pt_r[64]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[576] + u.imag[3][3] * amp_pt_r[576]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[64] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[576] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[64] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[576] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_8_7(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (256 << 1))
    {
        for (size_t tt = 0; tt < 256; tt += (128 << 1))
        {
            for (size_t ttt = 0; ttt < 128; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[128]   - u.imag[0][1] * amp_pt_i[128]) +
                           (u.real[0][2] * amp_pt_r[256]   - u.imag[0][2] * amp_pt_i[256]) +
                           (u.real[0][3] * amp_pt_r[384] - u.imag[0][3] * amp_pt_i[384]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[128]   - u.imag[1][1] * amp_pt_i[128]) +
                           (u.real[1][2] * amp_pt_r[256]   - u.imag[1][2] * amp_pt_i[256]) +
                           (u.real[1][3] * amp_pt_r[384] - u.imag[1][3] * amp_pt_i[384]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[128]   - u.imag[2][1] * amp_pt_i[128]) +
                           (u.real[2][2] * amp_pt_r[256]   - u.imag[2][2] * amp_pt_i[256]) +
                           (u.real[2][3] * amp_pt_r[384] - u.imag[2][3] * amp_pt_i[384]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[128]   - u.imag[3][1] * amp_pt_i[128]) +
                           (u.real[3][2] * amp_pt_r[256]   - u.imag[3][2] * amp_pt_i[256]) +
                           (u.real[3][3] * amp_pt_r[384] - u.imag[3][3] * amp_pt_i[384]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[128]   + u.imag[0][1] * amp_pt_r[128]) +
                           (u.real[0][2] * amp_pt_i[256]   + u.imag[0][2] * amp_pt_r[256]) +
                           (u.real[0][3] * amp_pt_i[384] + u.imag[0][3] * amp_pt_r[384]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[128]   + u.imag[1][1] * amp_pt_r[128]) +
                           (u.real[1][2] * amp_pt_i[256]   + u.imag[1][2] * amp_pt_r[256]) +
                           (u.real[1][3] * amp_pt_i[384] + u.imag[1][3] * amp_pt_r[384]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[128]   + u.imag[2][1] * amp_pt_r[128]) +
                           (u.real[2][2] * amp_pt_i[256]   + u.imag[2][2] * amp_pt_r[256]) +
                           (u.real[2][3] * amp_pt_i[384] + u.imag[2][3] * amp_pt_r[384]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[128]   + u.imag[3][1] * amp_pt_r[128]) +
                           (u.real[3][2] * amp_pt_i[256]   + u.imag[3][2] * amp_pt_r[256]) +
                           (u.real[3][3] * amp_pt_i[384] + u.imag[3][3] * amp_pt_r[384]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[128] = amp_r[1];
                amp_pt_r[256] = amp_r[2];
                amp_pt_r[384] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[128] = amp_i[1];
                amp_pt_i[256] = amp_i[2];
                amp_pt_i[384] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_7(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (128 << 1))
        {
            for (size_t ttt = 0; ttt < 128; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[128]   - u.imag[0][1] * amp_pt_i[128]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[640] - u.imag[0][3] * amp_pt_i[640]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[128]   - u.imag[1][1] * amp_pt_i[128]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[640] - u.imag[1][3] * amp_pt_i[640]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[128]   - u.imag[2][1] * amp_pt_i[128]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[640] - u.imag[2][3] * amp_pt_i[640]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[128]   - u.imag[3][1] * amp_pt_i[128]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[640] - u.imag[3][3] * amp_pt_i[640]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[128]   + u.imag[0][1] * amp_pt_r[128]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[640] + u.imag[0][3] * amp_pt_r[640]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[128]   + u.imag[1][1] * amp_pt_r[128]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[640] + u.imag[1][3] * amp_pt_r[640]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[128]   + u.imag[2][1] * amp_pt_r[128]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[640] + u.imag[2][3] * amp_pt_r[640]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[128]   + u.imag[3][1] * amp_pt_r[128]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[640] + u.imag[3][3] * amp_pt_r[640]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[128] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[640] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[128] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[640] = amp_i[3];
            } 
        }
    }
}


void two_qubit_gate_9_8(Statevector *sv, ComplexMatrix4 u) {
    double *amp_pt_r, *amp_pt_i;
    double amp_r[4], amp_i[4];

    for (size_t t = 0; t < sv->namp; t += (512 << 1))
    {
        for (size_t tt = 0; tt < 512; tt += (256 << 1))
        {
            for (size_t ttt = 0; ttt < 256; ttt++)
            {
                amp_pt_r = sv->real + t + tt + ttt;
                amp_pt_i = sv->imag + t + tt + ttt;
                amp_r[0] = (u.real[0][0] * amp_pt_r[0]   - u.imag[0][0] * amp_pt_i[0]) + 
                           (u.real[0][1] * amp_pt_r[256]   - u.imag[0][1] * amp_pt_i[256]) +
                           (u.real[0][2] * amp_pt_r[512]   - u.imag[0][2] * amp_pt_i[512]) +
                           (u.real[0][3] * amp_pt_r[768] - u.imag[0][3] * amp_pt_i[768]);

                amp_r[1] = (u.real[1][0] * amp_pt_r[0]   - u.imag[1][0] * amp_pt_i[0]) + 
                           (u.real[1][1] * amp_pt_r[256]   - u.imag[1][1] * amp_pt_i[256]) +
                           (u.real[1][2] * amp_pt_r[512]   - u.imag[1][2] * amp_pt_i[512]) +
                           (u.real[1][3] * amp_pt_r[768] - u.imag[1][3] * amp_pt_i[768]);

                amp_r[2] = (u.real[2][0] * amp_pt_r[0]   - u.imag[2][0] * amp_pt_i[0]) + 
                           (u.real[2][1] * amp_pt_r[256]   - u.imag[2][1] * amp_pt_i[256]) +
                           (u.real[2][2] * amp_pt_r[512]   - u.imag[2][2] * amp_pt_i[512]) +
                           (u.real[2][3] * amp_pt_r[768] - u.imag[2][3] * amp_pt_i[768]);

                amp_r[3] = (u.real[3][0] * amp_pt_r[0]   - u.imag[3][0] * amp_pt_i[0]) + 
                           (u.real[3][1] * amp_pt_r[256]   - u.imag[3][1] * amp_pt_i[256]) +
                           (u.real[3][2] * amp_pt_r[512]   - u.imag[3][2] * amp_pt_i[512]) +
                           (u.real[3][3] * amp_pt_r[768] - u.imag[3][3] * amp_pt_i[768]);

                amp_i[0] = (u.real[0][0] * amp_pt_i[0]   + u.imag[0][0] * amp_pt_r[0]) + 
                           (u.real[0][1] * amp_pt_i[256]   + u.imag[0][1] * amp_pt_r[256]) +
                           (u.real[0][2] * amp_pt_i[512]   + u.imag[0][2] * amp_pt_r[512]) +
                           (u.real[0][3] * amp_pt_i[768] + u.imag[0][3] * amp_pt_r[768]);

                amp_i[1] = (u.real[1][0] * amp_pt_i[0]   + u.imag[1][0] * amp_pt_r[0]) + 
                           (u.real[1][1] * amp_pt_i[256]   + u.imag[1][1] * amp_pt_r[256]) +
                           (u.real[1][2] * amp_pt_i[512]   + u.imag[1][2] * amp_pt_r[512]) +
                           (u.real[1][3] * amp_pt_i[768] + u.imag[1][3] * amp_pt_r[768]);

                amp_i[2] = (u.real[2][0] * amp_pt_i[0]   + u.imag[2][0] * amp_pt_r[0]) + 
                           (u.real[2][1] * amp_pt_i[256]   + u.imag[2][1] * amp_pt_r[256]) +
                           (u.real[2][2] * amp_pt_i[512]   + u.imag[2][2] * amp_pt_r[512]) +
                           (u.real[2][3] * amp_pt_i[768] + u.imag[2][3] * amp_pt_r[768]);

                amp_i[3] = (u.real[3][0] * amp_pt_i[0]   + u.imag[3][0] * amp_pt_r[0]) + 
                           (u.real[3][1] * amp_pt_i[256]   + u.imag[3][1] * amp_pt_r[256]) +
                           (u.real[3][2] * amp_pt_i[512]   + u.imag[3][2] * amp_pt_r[512]) +
                           (u.real[3][3] * amp_pt_i[768] + u.imag[3][3] * amp_pt_r[768]);

                amp_pt_r[0] = amp_r[0];
                amp_pt_r[256] = amp_r[1];
                amp_pt_r[512] = amp_r[2];
                amp_pt_r[768] = amp_r[3];
                amp_pt_i[0] = amp_i[0];
                amp_pt_i[256] = amp_i[1];
                amp_pt_i[512] = amp_i[2];
                amp_pt_i[768] = amp_i[3];
            } 
        }
    }
}


