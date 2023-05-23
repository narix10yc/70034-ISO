#include <stdio.h>
#include <stdlib.h>

// #include <smmintrin.h>

#include "timing.h"
#include "simd.h"


size_t N = 1 << 15;
double *x, *y, *rst;

void pre() {
    x = (double *) aligned_alloc(64, N * sizeof(double));
    y = (double *) aligned_alloc(64, N * sizeof(double));
    rst = (double *) aligned_alloc(64, N * sizeof(double));

    // printf("Memory allocation: \nx %p, y %p, rst %p\n", x, y, rst);
}

void post() {
    free(x);
    free(y);
    free(rst);
}

void no_sse_scalar() {
    no_sse_scalar_add_d(x, y, rst, N);
}

void scalar() {
    scalar_add_d(x, y, rst, N);
}

void sse() {
    sse_add_d(x, y, rst, N);
}

void avx() {
    avx_add_d(x, y, rst, N);
}

void avx512() {
    avx512_add_d(x, y, rst, N);
}


int main(int argc, char *argv[]) {
    // if (argc > 1) {
    //     if (atoi(argv[1]) > 25) {
    //         printf("Too large. Abort\n");
    //         return 0;
    //     }
    //     N = 1 << atoi(argv[1]);
    // }
    
    TimingResult tr;

    printf("No SSE scalar add\n");
    for (int n = 10; n < 25; n++)
    {
        N = 1 << n;
        timeit(pre, no_sse_scalar, post, &tr);
        printf("%.5g, ", tr.wall_stat[2] / tr.repetition * 1e6);
    }

    printf("\n\nScalar add\n");
    for (int n = 10; n < 25; n++)
    {
        N = 1 << n;
        timeit(pre, scalar, post, &tr);
        printf("%.5g, ", tr.wall_stat[2] / tr.repetition * 1e6);
    }
    
    printf("\n\nSSE add\n");
    for (int n = 10; n < 25; n++)
    {
        N = 1 << n;
        timeit(pre, sse, post, &tr);
        printf("%.5g, ", tr.wall_stat[2] / tr.repetition * 1e6);
    }
    
    printf("\n\nAVX add\n");
    for (int n = 10; n < 25; n++)
    {
        N = 1 << n;
        timeit(pre, avx, post, &tr);
        printf("%.5g, ", tr.wall_stat[2] / tr.repetition * 1e6);
    }
    
    printf("\n\nAVX512 add\n");
    for (int n = 10; n < 25; n++)
    {
        N = 1 << n;
        timeit(pre, avx512, post, &tr);
        printf("%.5g, ", tr.wall_stat[2] / tr.repetition * 1e6);
    }
    
    printf("\n");


    return 0;
}

