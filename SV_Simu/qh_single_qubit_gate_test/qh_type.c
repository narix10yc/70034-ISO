// #include <stdio.h>
#include <stdlib.h>
#include "qh_type.h"

Camp *create_amp(int n) {
    size_t N = 1 << n;
    double *real = (double *) aligned_alloc(64, N * sizeof(double));
    double *imag = (double *) aligned_alloc(64, N * sizeof(double));

    Camp *amp = (Camp *) malloc(sizeof(Camp));
    amp->real = real;
    amp->imag = imag;

    return amp;
}

void destroy_amp(Camp *amp) {
    free(amp->real);
    free(amp->imag);
    free(amp);
}
