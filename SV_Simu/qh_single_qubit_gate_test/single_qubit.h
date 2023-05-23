#ifndef SINGLE_QUBIT_H
#define SINGLE_QUBIT_H

#include <stdlib.h>
#include "qh_type.h"


void extract_indices(int n, int k, size_t *a, size_t *b);

void apply_on_k(int n, int k, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);

void apply_on_k_prompt(int n, int k, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void apply_on_k_prompt_simd(int n, int k, Camp amp, ComplexMatrix2 u);
void apply_on_k_prompt_matrix(int n, int k, Camp amp, ComplexMatrix2 u);
void apply_on_k_prompt_multi_threading(int n, int k, Camp amp, ComplexMatrix2 u);


void explicit_on_0(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_1(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_2(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_3(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_4(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_5(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_6(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_7(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_8(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_9(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);
void explicit_on_10(int n, Camp amp, Cdouble a, Cdouble b, Cdouble c, Cdouble d);

void apply_on_1_manual(int n, double *amp, double a, double b, double c);


void apply_on_2_manual(int n, double *amp, double a, double b, double c);


void apply_on_3_manual(int n, double *amp, double a, double b, double c);


void apply_on_4_manual(int n, double *amp, double a, double b, double c);

void statevec_unitaryLocal(Qureg qureg, int targetQubit, ComplexMatrix2 u);

#endif