#include <stdlib.h>
#include <stdio.h>

#include "timing.h"
#include "complex_mul.h"

double *xr, *xi, *yr, *yi, *rstr, *rsti;
double *x, *y, *rst;
Cdouble *xc, *yc, *rstc;
complex double *xch, *ych, *rstch;

int n = 18;
size_t N;

void pre_separate() {
    N = 1 << n;
    xr = (double *) aligned_alloc(64, N * sizeof(double));
    xi = (double *) aligned_alloc(64, N * sizeof(double));
    yr = (double *) aligned_alloc(64, N * sizeof(double));
    yi = (double *) aligned_alloc(64, N * sizeof(double));
    rstr = (double *) aligned_alloc(64, N * sizeof(double));
    rsti = (double *) aligned_alloc(64, N * sizeof(double));
}

void post_separate() {
    free(xr);
    free(xi);
    free(yr);
    free(yi);
    free(rstr);
    free(rsti);
}

void pre_combined() {
    N = 1 << n;
    x = (double *) aligned_alloc(64, 2*N * sizeof(double));
    y = (double *) aligned_alloc(64, 2*N * sizeof(double));
    rst = (double *) aligned_alloc(64, 2*N * sizeof(double));
}

void post_combined() {
    free(x);
    free(y);
    free(rst);
}

void pre_cdouble() {
    N = 1 << n;
    xc = (Cdouble *) aligned_alloc(64, N * sizeof(Cdouble));
    yc = (Cdouble *) aligned_alloc(64, N * sizeof(Cdouble));
    rstc = (Cdouble *) aligned_alloc(64, N * sizeof(Cdouble));
}

void post_cdouble() {
    free(xc);
    free(yc);
    free(rstc);
}

void pre_complexh() {
    N = 1 << n;
    xch = (complex double *) aligned_alloc(64, N * sizeof(complex double));
    ych = (complex double *) aligned_alloc(64, N * sizeof(complex double));
    rstch = (complex double *) aligned_alloc(64, N * sizeof(complex double));
}

void post_complexh() {
    free(xch);
    free(ych);
    free(rstch);
}

void timeit_scalar_separate() {
    scalar_separate(xr, xi, yr, yi, rstr, rsti, N);
}

void timeit_scalar_combined() {
    scalar_combined(x, y, rst, N);
}

void timeit_scalar_cdouble() {
    scalar_cdouble(xc, yc, rstc, N);
}

void timeit_scalar_complexh() {
    scalar_complexh(xch, ych, rstch, N);
}


void timeit_simd_separate() {
    simd_separate(xr, xi, yr, yi, rstr, rsti, N);
}

void timeit_simd_combined() {
    simd_combined(x, y, rst, N);
}

void timeit_simd_cdouble() {
    simd_cdouble(xc, yc, rstc, N);
}


void main(int argc, char *argv[]) {
    if (argc > 1) {
        n = atoi(argv[1]);
    }

    TimingResult tr;

    printf("Scalar Separate:\n");
    timeit(pre_separate, timeit_scalar_separate, post_separate, &tr);
    print_timing_result(&tr);


    printf("\n\nScalar Combined:\n");
    timeit(pre_combined, timeit_scalar_combined, post_combined, &tr);
    print_timing_result(&tr);


    printf("\n\nScalar Cdouble:\n");
    timeit(pre_cdouble, timeit_scalar_cdouble, post_cdouble, &tr);
    print_timing_result(&tr);


    printf("\n\nScalar complex.h:\n");
    timeit(pre_complexh, timeit_scalar_complexh, post_complexh, &tr);
    print_timing_result(&tr);


    // printf("\n\nSIMD Separate:\n");
    // timeit(pre_separate, timeit_simd_separate, post_separate, &tr);
    // print_timing_result(&tr);


    // printf("\n\nSIMD Combined:\n");
    // timeit(pre_combined, timeit_simd_combined, post_combined, &tr);
    // print_timing_result(&tr);


    // printf("\n\nSIMD Cdouble:\n");
    // timeit(pre_cdouble, timeit_simd_cdouble, post_cdouble, &tr);
    // print_timing_result(&tr);

}