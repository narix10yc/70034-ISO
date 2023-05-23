#include "timing.h"
#include "two_qubit_manual.h"
#include "qh_type.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #include <omp.h>

size_t n = 19;
Statevector *sv;
ComplexMatrix4 gate = (ComplexMatrix4) {{{0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}, 
                                         {0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}, 
                                         {0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}}};

void pre() {
    size_t N = 1 << n;
    sv = svCreate(n);
    svRandomize(sv);
    
}

void post() {
    svDestroy(sv);
}



int main(int argc, char *argv[]) {

    if (argc > 1) {
        n = atoi(argv[1]);
    }

    TimingResult rst;

    // printf("General:\n");
    // for (int k = 0; k < 10; k++)
    // {
    //     timeit(pre, timeit_general[k], post, &rst);
    //     printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    // }


    printf("\n");

}

