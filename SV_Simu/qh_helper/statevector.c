#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "statevector.h"

Statevector *svCreate(int nqubits) {
    size_t namp = 1 << nqubits;
    double *real = (double *) aligned_alloc(64, namp * sizeof(double));
    double *imag = (double *) aligned_alloc(64, namp * sizeof(double));
    Statevector *sv = (Statevector *) malloc(sizeof(Statevector));
    *sv = (Statevector) { .nqubits = nqubits, .namp = namp, .real = real, .imag = imag };
    return sv;
}

void svDestroy(Statevector *sv) {
    free(sv->real);
    free(sv->imag);
    free(sv);
}

void svInit(Statevector *sv) {
    for (size_t i = 0; i < sv->namp; i++) 
    {
        sv->real[i] = 0;
        sv->imag[i] = 0;
    }
    sv->real[0] = 1;
}

void svRandomize(Statevector *sv) {
    srand(time(NULL));
    double sum = 0;
    for (size_t i = 0; i < sv->namp; i++) 
    {
        sv->real[i] = (double)rand() / RAND_MAX;
        sv->imag[i] = (double)rand() / RAND_MAX;

        sum += (sv->real[i] * sv->real[i]) + (sv->imag[i] * sv->imag[i]);
    }
    
    for (size_t i = 0; i < sv->namp; i++) 
    {
        sv->real[i] /= sum;
        sv->imag[i] /= sum;
    }
}

void svPrintAmplitude(Statevector *sv) {
    for (size_t i = 0; i < sv->namp; i += 4)
    {
        for (size_t ii = 0; ii < 4; ii++)
        {
            printf("%.4f + %.4fi  ", sv->real[i+ii], sv->imag[i+ii]);
        }
        printf("\n");
    }
}

void svPrintProbability(Statevector *sv) {
    double prob;
    for (size_t i = 0; i < sv->namp; i += 8)
    {
        for (size_t ii = 0; ii < 8; ii++)
        {
            prob = (sv->real[i+ii] * sv->real[i+ii]) + (sv->imag[i+ii] * sv->imag[i+ii]);
            printf("%.4f ", prob);
        }
        printf("\n");
    }
}
