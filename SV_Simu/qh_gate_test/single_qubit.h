#ifndef SINGLE_QUBIT_H
#define SINGLE_QUBIT_H

#include <stdlib.h>

void extract_indices(int n, int k, size_t *a, size_t *b);

void apply_on_k(int n, int k, double *amp, double a, double b, double c);

void apply_on_k_prompt(int n, int k, double *amp, double a, double b, double c);

void apply_on_0(int n, double *amp, double a, double b, double c);

void apply_on_1(int n, double *amp, double a, double b, double c);

void apply_on_1_manual(int n, double *amp, double a, double b, double c);

void apply_on_2(int n, double *amp, double a, double b, double c);

void apply_on_2_manual(int n, double *amp, double a, double b, double c);

void apply_on_3(int n, double *amp, double a, double b, double c);

void apply_on_3_manual(int n, double *amp, double a, double b, double c);

void apply_on_4(int n, double *amp, double a, double b, double c);

void apply_on_4_manual(int n, double *amp, double a, double b, double c);

#endif