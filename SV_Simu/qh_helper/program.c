#include <stdlib.h>
#include <stdio.h>

#include "timing.h"
#include "qh_gates.h"
#include "statevector.h"

Statevector *sv;
int n = 19;
int k = 0;
ComplexMatrix2 gate = (ComplexMatrix2) {{{0.5, 0.5}, {0.5, 0.5}}, {{0.5, 0.5}, {0.5, 0.5}}};

void pre() {
    sv = svCreate(n);
    svRandomize(sv);
}

void post() {
    svDestroy(sv);
}

void timeit_gate() {
    applyGateSingleUnitary(sv, k, gate);
}

int main(int argc, char *argv[]) {
    if (argc > 1) {
        k = atoi(argv[1]);
    }

    TimingResult tr;
    timeit(pre, timeit_gate, post, &tr);
    print_timing_result(&tr);


    return 0;
}