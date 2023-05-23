#ifndef QH_TYPE_H
#define QH_TYPE_H

typedef struct Cdouble {
    double real, imag;
} Cdouble;


typedef struct Camp {
    double *real;
    double *imag;
} Camp;

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

typedef struct Qureg
{
    long long int numAmpsPerChunk;
    double *real; 
    double *imag;
} Qureg;

Camp *create_amp(int n);

void destroy_amp(Camp *amp);


#endif