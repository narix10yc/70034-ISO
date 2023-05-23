#include "timing.h"
#include "single_qubit.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

size_t n = 20;
double *amp;

void pre() {
    amp = (double *) aligned_alloc(64, (1 << n) * sizeof(double));
    srand(time(NULL));
    double sum = 0;
    for (size_t i = 0; i < (1 << n); i++)
    {
        amp[i] = rand();
        sum += amp[i] * amp[i];
    }
    
    for (size_t i = 0; i < (1 << n); i++)
    {
        amp[i] /= sum;
    }
    
}

void post() {
    free(amp);
}

// Qubit 0
void timeit_apply_on_k0() {
    apply_on_k(n, 0, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_k0_prompt() {
    apply_on_k_prompt(n, 0, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_0() {
    apply_on_0(n, amp, 0.5, 0.5, 0.5);
}

// Qubit 1
void timeit_apply_on_k1() {
    apply_on_k(n, 1, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_k1_prompt() {
    apply_on_k_prompt(n, 1, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_1() {
    apply_on_1(n, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_1_manual() {
    apply_on_1_manual(n, amp, 0.5, 0.5, 0.5);
}

// Qubit 2
void timeit_apply_on_k2() {
    apply_on_k(n, 2, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_k2_prompt() {
    apply_on_k_prompt(n, 2, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_2() {
    apply_on_2(n, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_2_manual() {
    apply_on_2_manual(n, amp, 0.5, 0.5, 0.5);
}

// Qubit 3
void timeit_apply_on_k3() {
    apply_on_k(n, 3, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_k3_prompt() {
    apply_on_k_prompt(n, 3, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_3() {
    apply_on_3(n, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_3_manual() {
    apply_on_3_manual(n, amp, 0.5, 0.5, 0.5);
}

// Qubit 3
void timeit_apply_on_k4() {
    apply_on_k(n, 4, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_k4_prompt() {
    apply_on_k_prompt(n, 4, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_4() {
    apply_on_4(n, amp, 0.5, 0.5, 0.5);
}

void timeit_apply_on_4_manual() {
    apply_on_4_manual(n, amp, 0.5, 0.5, 0.5);
}


int main(int argc, char *argv[]) {

    if (argc > 1) {
        n = atoi(argv[1]);
    }

    TimingResult rst;

    printf("Apply on 0:\n");
    timeit(pre, timeit_apply_on_k0, post, &rst);
    printf("\tGeneral: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_k0_prompt, post, &rst);
    printf("\tPrompt: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_0, post, &rst);
    printf("\tSpecial: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);

    printf("Apply on 1:\n");
    timeit(pre, timeit_apply_on_k1, post, &rst);
    printf("\tGeneral: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_k1_prompt, post, &rst);
    printf("\tPrompt: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_1, post, &rst);
    printf("\tSpecial: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_1_manual, post, &rst);
    printf("\tManual: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);

    printf("Apply on 2:\n");
    timeit(pre, timeit_apply_on_k2, post, &rst);
    printf("\tGeneral: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_k2_prompt, post, &rst);
    printf("\tPrompt: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_2, post, &rst);
    printf("\tSpecial: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_2_manual, post, &rst);
    printf("\tManual: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    
    printf("Apply on 3:\n");
    timeit(pre, timeit_apply_on_k3, post, &rst);
    printf("\tGeneral: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_k3_prompt, post, &rst);
    printf("\tPrompt: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_3, post, &rst);
    printf("\tSpecial: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_3_manual, post, &rst);
    printf("\tManual: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);

    printf("Apply on 4:\n");
    timeit(pre, timeit_apply_on_k4, post, &rst);
    printf("\tGeneral: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_k4_prompt, post, &rst);
    printf("\tPrompt: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_4, post, &rst);
    printf("\tSpecial: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);
    timeit(pre, timeit_apply_on_4_manual, post, &rst);
    printf("\tManual: %f\n", rst.cpu_stat[2] / rst.repetition * 1e6);

}

