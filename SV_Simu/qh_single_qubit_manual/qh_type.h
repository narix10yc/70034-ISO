#ifndef QH_TYPE_H
#define QH_TYPE_H

#include <stdlib.h>

typedef struct Statevector {
    int nqubits;
    size_t namp;
    double *real;
    double *imag;
} Statevector;


typedef struct ComplexMatrix2
{
    double real[2][2];
    double imag[2][2];
} ComplexMatrix2;


typedef struct ComplexMatrix4
{
    double real[4][4];
    double imag[4][4];
} ComplexMatrix4;


#endif