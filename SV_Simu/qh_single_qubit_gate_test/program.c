#include "timing.h"
#include "single_qubit.h"
#include "qh_type.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <omp.h>

size_t n = 16;
Camp *amp;
Qureg qureg;
ComplexMatrix2 gate = (ComplexMatrix2) {{{0.5, 0.5}, {0.5, 0.5}}, {{0.5, 0.5}, {0.5, 0.5}}};


void pre() {
    size_t N = 1 << n;
    amp = create_amp(n);
    srand(time(NULL));
    double sum = 0;
    for (size_t i = 0; i < (1 << n); i++)
    {
        amp->real[i] = rand();
        amp->imag[i] = rand();

        sum += (amp->real[i] * amp->real[i]) + (amp->imag[i] * amp->imag[i]);
    }
    
    for (size_t i = 0; i < (1 << n); i++)
    {
        amp->real[i] /= sum;
        amp->imag[i] /= sum;
    }

    double *real = (double *) aligned_alloc(64, N * sizeof(double));
    double *imag = (double *) aligned_alloc(64, N * sizeof(double));

    qureg = (Qureg) {1 << n, real, imag};
    
}

void post() {
    destroy_amp(amp);
    free(qureg.real);
    free(qureg.imag);
}


// Explicit Method
void timeit_explicit_on_0() {
    explicit_on_0(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_1() {
    explicit_on_1(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_2() {
    explicit_on_2(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_3() {
    explicit_on_3(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_4() {
    explicit_on_4(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_5() {
    explicit_on_5(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_6() {
    explicit_on_6(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_7() {
    explicit_on_7(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_8() {
    explicit_on_8(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_explicit_on_9() {
    explicit_on_9(n, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}


void (*timeit_explicit[10])() = 
{ timeit_explicit_on_0, timeit_explicit_on_1, timeit_explicit_on_2, timeit_explicit_on_3, timeit_explicit_on_4,
  timeit_explicit_on_5, timeit_explicit_on_6, timeit_explicit_on_7, timeit_explicit_on_8, timeit_explicit_on_9 };


// General Method
void timeit_general_on_0() {
    apply_on_k(n, 0, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_1() {
    apply_on_k(n, 1, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_2() {
    apply_on_k(n, 2, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_3() {
    apply_on_k(n, 3, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_4() {
    apply_on_k(n, 4, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_5() {
    apply_on_k(n, 5, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_6() {
    apply_on_k(n, 6, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_7() {
    apply_on_k(n, 7, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_8() {
    apply_on_k(n, 8, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}
void timeit_general_on_9() {
    apply_on_k(n, 9, *amp, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5}, (Cdouble){0.5, 0.5});
}


void (*timeit_general[10])() = 
{ timeit_general_on_0, timeit_general_on_1, timeit_general_on_2, timeit_general_on_3, timeit_general_on_4,
  timeit_general_on_5, timeit_general_on_6, timeit_general_on_7, timeit_general_on_8, timeit_general_on_9 };


// Prompt Method
void timeit_prompt_on_0() {
    apply_on_k_prompt_multi_threading(n, 0, *amp, gate);
}
void timeit_prompt_on_1() {
    apply_on_k_prompt_multi_threading(n, 1, *amp, gate);
}
void timeit_prompt_on_2() {
    apply_on_k_prompt_multi_threading(n, 2, *amp, gate);
}
void timeit_prompt_on_3() {
    apply_on_k_prompt_multi_threading(n, 3, *amp, gate);
}
void timeit_prompt_on_4() {
    apply_on_k_prompt_multi_threading(n, 4, *amp, gate);
}
void timeit_prompt_on_5() {
    apply_on_k_prompt_multi_threading(n, 5, *amp, gate);
}
void timeit_prompt_on_6() {
    apply_on_k_prompt_multi_threading(n, 6, *amp, gate);
}
void timeit_prompt_on_7() {
    apply_on_k_prompt_multi_threading(n, 7, *amp, gate);
}
void timeit_prompt_on_8() {
    apply_on_k_prompt_multi_threading(n, 8, *amp, gate);
}
void timeit_prompt_on_9() {
    apply_on_k_prompt_multi_threading(n, 9, *amp, gate);
}
void timeit_prompt_on_10() {
    apply_on_k_prompt_multi_threading(n, 10, *amp, gate);
}
void timeit_prompt_on_11() {
    apply_on_k_prompt_multi_threading(n, 11, *amp, gate);
}
void timeit_prompt_on_12() {
    apply_on_k_prompt_multi_threading(n, 12, *amp, gate);
}
void timeit_prompt_on_13() {
    apply_on_k_prompt_multi_threading(n, 13, *amp, gate);
}
void timeit_prompt_on_14() {
    apply_on_k_prompt_multi_threading(n, 14, *amp, gate);
}
void timeit_prompt_on_15() {
    apply_on_k_prompt_multi_threading(n, 15, *amp, gate);
}
void timeit_prompt_on_16() {
    apply_on_k_prompt_multi_threading(n, 16, *amp, gate);
}
void timeit_prompt_on_17() {
    apply_on_k_prompt_multi_threading(n, 17, *amp, gate);
}


void (*timeit_prompt[18])() = 
{ timeit_prompt_on_0, timeit_prompt_on_1, timeit_prompt_on_2, timeit_prompt_on_3, timeit_prompt_on_4,
  timeit_prompt_on_5, timeit_prompt_on_6, timeit_prompt_on_7, timeit_prompt_on_8, timeit_prompt_on_9,
  timeit_prompt_on_10, timeit_prompt_on_11, timeit_prompt_on_12, timeit_prompt_on_13, timeit_prompt_on_14,
  timeit_prompt_on_15, timeit_prompt_on_16, timeit_prompt_on_17
};



// QuEST 
void timeit_QuEST_0() {
    statevec_unitaryLocal(qureg, 0, gate);
}
void timeit_QuEST_1() {
    statevec_unitaryLocal(qureg, 1, gate);
}
void timeit_QuEST_2() {
    statevec_unitaryLocal(qureg, 2, gate);
}
void timeit_QuEST_3() {
    statevec_unitaryLocal(qureg, 3, gate);
}
void timeit_QuEST_4() {
    statevec_unitaryLocal(qureg, 4, gate);
}
void timeit_QuEST_5() {
    statevec_unitaryLocal(qureg, 5, gate);
}
void timeit_QuEST_6() {
    statevec_unitaryLocal(qureg, 6, gate);
}
void timeit_QuEST_7() {
    statevec_unitaryLocal(qureg, 7, gate);
}
void timeit_QuEST_8() {
    statevec_unitaryLocal(qureg, 8, gate);
}
void timeit_QuEST_9() {
    statevec_unitaryLocal(qureg, 9, gate);
}
void timeit_QuEST_10() {
    statevec_unitaryLocal(qureg, 10, gate);
}
void timeit_QuEST_11() {
    statevec_unitaryLocal(qureg, 11, gate);
}
void timeit_QuEST_12() {
    statevec_unitaryLocal(qureg, 12, gate);
}
void timeit_QuEST_13() {
    statevec_unitaryLocal(qureg, 13, gate);
}
void timeit_QuEST_14() {
    statevec_unitaryLocal(qureg, 14, gate);
}
void timeit_QuEST_15() {
    statevec_unitaryLocal(qureg, 15, gate);
}
void timeit_QuEST_16() {
    statevec_unitaryLocal(qureg, 16, gate);
}
void timeit_QuEST_17() {
    statevec_unitaryLocal(qureg, 17, gate);
}


void (*timeit_QuEST[18])() = 
{ timeit_QuEST_0, timeit_QuEST_1, timeit_QuEST_2, timeit_QuEST_3, timeit_QuEST_4,
  timeit_QuEST_5, timeit_QuEST_6, timeit_QuEST_7, timeit_QuEST_8, timeit_QuEST_9,
  timeit_QuEST_10, timeit_QuEST_11, timeit_QuEST_12, timeit_QuEST_13, timeit_QuEST_14,
  timeit_QuEST_15, timeit_QuEST_16, timeit_QuEST_17
};


int main(int argc, char *argv[]) {

    int num_threads = 1;
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }

    omp_set_num_threads(num_threads);

    TimingResult rst;

    // printf("General:\n");
    // for (int k = 0; k < 10; k++)
    // {
    //     timeit(pre, timeit_general[k], post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    // printf("\n\nPrompt:\n");
    // for (int k = 0; k < 18; k++)
    // {
    //     timeit(pre, timeit_prompt[k], post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }

    // printf("\n\nExplicit:\n");
    // for (int k = 0; k < 10; k++)
    // {
    //     timeit(pre, timeit_explicit[k], post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }
    

    // printf("\n\nQuEST:\n");
    // for (int k = 0; k < 18; k++)
    // {
    //     timeit(pre, timeit_QuEST[k], post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }


    int num_threads_list[] = {1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32};
    printf("\n\nPrompt Multi-Threading:\n");
    for (int nt = 0; nt < 13; nt++)
    {
        sleep(1);
        omp_set_num_threads(num_threads_list[nt]);
        timeit(pre, timeit_prompt[n-1], post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    }  

    printf("\n\nQuEST Multi-Threading:\n");
    for (int nt = 0; nt < 13; nt++)
    {
        sleep(1);
        omp_set_num_threads(num_threads_list[nt]);
        timeit(pre, timeit_QuEST[n-1], post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    }  
    
    printf("\n");

}

