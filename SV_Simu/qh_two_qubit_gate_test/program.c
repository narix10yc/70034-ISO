#include "timing.h"

#include "two_qubit_explicit.h"
#include "two_qubit_general.h"
#include "two_qubit_quest.h"

#include "statevector.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <unistd.h>

#include <omp.h>

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

void timeit_general_1_0() { general_two_qubit_gate(sv, 1, 0, gate); }
void timeit_general_2_0() { general_two_qubit_gate(sv, 2, 0, gate); }
void timeit_general_3_0() { general_two_qubit_gate(sv, 3, 0, gate); }
void timeit_general_4_0() { general_two_qubit_gate(sv, 4, 0, gate); }
void timeit_general_5_0() { general_two_qubit_gate(sv, 5, 0, gate); }
void timeit_general_6_0() { general_two_qubit_gate(sv, 6, 0, gate); }
void timeit_general_7_0() { general_two_qubit_gate(sv, 7, 0, gate); }
void timeit_general_8_0() { general_two_qubit_gate(sv, 8, 0, gate); }
void timeit_general_9_0() { general_two_qubit_gate(sv, 9, 0, gate); }
void timeit_general_2_1() { general_two_qubit_gate(sv, 2, 1, gate); }
void timeit_general_3_1() { general_two_qubit_gate(sv, 3, 1, gate); }
void timeit_general_4_1() { general_two_qubit_gate(sv, 4, 1, gate); }
void timeit_general_5_1() { general_two_qubit_gate(sv, 5, 1, gate); }
void timeit_general_6_1() { general_two_qubit_gate(sv, 6, 1, gate); }
void timeit_general_7_1() { general_two_qubit_gate(sv, 7, 1, gate); }
void timeit_general_8_1() { general_two_qubit_gate(sv, 8, 1, gate); }
void timeit_general_9_1() { general_two_qubit_gate(sv, 9, 1, gate); }
void timeit_general_3_2() { general_two_qubit_gate(sv, 3, 2, gate); }
void timeit_general_4_2() { general_two_qubit_gate(sv, 4, 2, gate); }
void timeit_general_5_2() { general_two_qubit_gate(sv, 5, 2, gate); }
void timeit_general_6_2() { general_two_qubit_gate(sv, 6, 2, gate); }
void timeit_general_7_2() { general_two_qubit_gate(sv, 7, 2, gate); }
void timeit_general_8_2() { general_two_qubit_gate(sv, 8, 2, gate); }
void timeit_general_9_2() { general_two_qubit_gate(sv, 9, 2, gate); }
void timeit_general_4_3() { general_two_qubit_gate(sv, 4, 3, gate); }
void timeit_general_5_3() { general_two_qubit_gate(sv, 5, 3, gate); }
void timeit_general_6_3() { general_two_qubit_gate(sv, 6, 3, gate); }
void timeit_general_7_3() { general_two_qubit_gate(sv, 7, 3, gate); }
void timeit_general_8_3() { general_two_qubit_gate(sv, 8, 3, gate); }
void timeit_general_9_3() { general_two_qubit_gate(sv, 9, 3, gate); }
void timeit_general_5_4() { general_two_qubit_gate(sv, 5, 4, gate); }
void timeit_general_6_4() { general_two_qubit_gate(sv, 6, 4, gate); }
void timeit_general_7_4() { general_two_qubit_gate(sv, 7, 4, gate); }
void timeit_general_8_4() { general_two_qubit_gate(sv, 8, 4, gate); }
void timeit_general_9_4() { general_two_qubit_gate(sv, 9, 4, gate); }
void timeit_general_6_5() { general_two_qubit_gate(sv, 6, 5, gate); }
void timeit_general_7_5() { general_two_qubit_gate(sv, 7, 5, gate); }
void timeit_general_8_5() { general_two_qubit_gate(sv, 8, 5, gate); }
void timeit_general_9_5() { general_two_qubit_gate(sv, 9, 5, gate); }
void timeit_general_7_6() { general_two_qubit_gate(sv, 7, 6, gate); }
void timeit_general_8_6() { general_two_qubit_gate(sv, 8, 6, gate); }
void timeit_general_9_6() { general_two_qubit_gate(sv, 9, 6, gate); }
void timeit_general_8_7() { general_two_qubit_gate(sv, 8, 7, gate); }
void timeit_general_9_7() { general_two_qubit_gate(sv, 9, 7, gate); }
void timeit_general_9_8() { general_two_qubit_gate(sv, 9, 8, gate); }


void (*timeit_general[45])() = {
	timeit_general_1_0,
	timeit_general_2_0,
	timeit_general_3_0,
	timeit_general_4_0,
	timeit_general_5_0,
	timeit_general_6_0,
	timeit_general_7_0,
	timeit_general_8_0,
	timeit_general_9_0,
	timeit_general_2_1,
	timeit_general_3_1,
	timeit_general_4_1,
	timeit_general_5_1,
	timeit_general_6_1,
	timeit_general_7_1,
	timeit_general_8_1,
	timeit_general_9_1,
	timeit_general_3_2,
	timeit_general_4_2,
	timeit_general_5_2,
	timeit_general_6_2,
	timeit_general_7_2,
	timeit_general_8_2,
	timeit_general_9_2,
	timeit_general_4_3,
	timeit_general_5_3,
	timeit_general_6_3,
	timeit_general_7_3,
	timeit_general_8_3,
	timeit_general_9_3,
	timeit_general_5_4,
	timeit_general_6_4,
	timeit_general_7_4,
	timeit_general_8_4,
	timeit_general_9_4,
	timeit_general_6_5,
	timeit_general_7_5,
	timeit_general_8_5,
	timeit_general_9_5,
	timeit_general_7_6,
	timeit_general_8_6,
	timeit_general_9_6,
	timeit_general_8_7,
	timeit_general_9_7,
	timeit_general_9_8,
};


void timeit_explicit_1_0() { two_qubit_gate_1_0(sv, gate); }
void timeit_explicit_2_0() { two_qubit_gate_2_0(sv, gate); }
void timeit_explicit_3_0() { two_qubit_gate_3_0(sv, gate); }
void timeit_explicit_4_0() { two_qubit_gate_4_0(sv, gate); }
void timeit_explicit_5_0() { two_qubit_gate_5_0(sv, gate); }
void timeit_explicit_6_0() { two_qubit_gate_6_0(sv, gate); }
void timeit_explicit_7_0() { two_qubit_gate_7_0(sv, gate); }
void timeit_explicit_8_0() { two_qubit_gate_8_0(sv, gate); }
void timeit_explicit_9_0() { two_qubit_gate_9_0(sv, gate); }
void timeit_explicit_2_1() { two_qubit_gate_2_1(sv, gate); }
void timeit_explicit_3_1() { two_qubit_gate_3_1(sv, gate); }
void timeit_explicit_4_1() { two_qubit_gate_4_1(sv, gate); }
void timeit_explicit_5_1() { two_qubit_gate_5_1(sv, gate); }
void timeit_explicit_6_1() { two_qubit_gate_6_1(sv, gate); }
void timeit_explicit_7_1() { two_qubit_gate_7_1(sv, gate); }
void timeit_explicit_8_1() { two_qubit_gate_8_1(sv, gate); }
void timeit_explicit_9_1() { two_qubit_gate_9_1(sv, gate); }
void timeit_explicit_3_2() { two_qubit_gate_3_2(sv, gate); }
void timeit_explicit_4_2() { two_qubit_gate_4_2(sv, gate); }
void timeit_explicit_5_2() { two_qubit_gate_5_2(sv, gate); }
void timeit_explicit_6_2() { two_qubit_gate_6_2(sv, gate); }
void timeit_explicit_7_2() { two_qubit_gate_7_2(sv, gate); }
void timeit_explicit_8_2() { two_qubit_gate_8_2(sv, gate); }
void timeit_explicit_9_2() { two_qubit_gate_9_2(sv, gate); }
void timeit_explicit_4_3() { two_qubit_gate_4_3(sv, gate); }
void timeit_explicit_5_3() { two_qubit_gate_5_3(sv, gate); }
void timeit_explicit_6_3() { two_qubit_gate_6_3(sv, gate); }
void timeit_explicit_7_3() { two_qubit_gate_7_3(sv, gate); }
void timeit_explicit_8_3() { two_qubit_gate_8_3(sv, gate); }
void timeit_explicit_9_3() { two_qubit_gate_9_3(sv, gate); }
void timeit_explicit_5_4() { two_qubit_gate_5_4(sv, gate); }
void timeit_explicit_6_4() { two_qubit_gate_6_4(sv, gate); }
void timeit_explicit_7_4() { two_qubit_gate_7_4(sv, gate); }
void timeit_explicit_8_4() { two_qubit_gate_8_4(sv, gate); }
void timeit_explicit_9_4() { two_qubit_gate_9_4(sv, gate); }
void timeit_explicit_6_5() { two_qubit_gate_6_5(sv, gate); }
void timeit_explicit_7_5() { two_qubit_gate_7_5(sv, gate); }
void timeit_explicit_8_5() { two_qubit_gate_8_5(sv, gate); }
void timeit_explicit_9_5() { two_qubit_gate_9_5(sv, gate); }
void timeit_explicit_7_6() { two_qubit_gate_7_6(sv, gate); }
void timeit_explicit_8_6() { two_qubit_gate_8_6(sv, gate); }
void timeit_explicit_9_6() { two_qubit_gate_9_6(sv, gate); }
void timeit_explicit_8_7() { two_qubit_gate_8_7(sv, gate); }
void timeit_explicit_9_7() { two_qubit_gate_9_7(sv, gate); }
void timeit_explicit_9_8() { two_qubit_gate_9_8(sv, gate); }


void (*timeit_explicit[45])() = {
	timeit_explicit_1_0,
	timeit_explicit_2_0,
	timeit_explicit_3_0,
	timeit_explicit_4_0,
	timeit_explicit_5_0,
	timeit_explicit_6_0,
	timeit_explicit_7_0,
	timeit_explicit_8_0,
	timeit_explicit_9_0,
	timeit_explicit_2_1,
	timeit_explicit_3_1,
	timeit_explicit_4_1,
	timeit_explicit_5_1,
	timeit_explicit_6_1,
	timeit_explicit_7_1,
	timeit_explicit_8_1,
	timeit_explicit_9_1,
	timeit_explicit_3_2,
	timeit_explicit_4_2,
	timeit_explicit_5_2,
	timeit_explicit_6_2,
	timeit_explicit_7_2,
	timeit_explicit_8_2,
	timeit_explicit_9_2,
	timeit_explicit_4_3,
	timeit_explicit_5_3,
	timeit_explicit_6_3,
	timeit_explicit_7_3,
	timeit_explicit_8_3,
	timeit_explicit_9_3,
	timeit_explicit_5_4,
	timeit_explicit_6_4,
	timeit_explicit_7_4,
	timeit_explicit_8_4,
	timeit_explicit_9_4,
	timeit_explicit_6_5,
	timeit_explicit_7_5,
	timeit_explicit_8_5,
	timeit_explicit_9_5,
	timeit_explicit_7_6,
	timeit_explicit_8_6,
	timeit_explicit_9_6,
	timeit_explicit_8_7,
	timeit_explicit_9_7,
	timeit_explicit_9_8,
};

void timeit_quest_1_0() { quest_two_qubit_gate(sv, 1, 0, gate); }
void timeit_quest_2_0() { quest_two_qubit_gate(sv, 2, 0, gate); }
void timeit_quest_3_0() { quest_two_qubit_gate(sv, 3, 0, gate); }
void timeit_quest_4_0() { quest_two_qubit_gate(sv, 4, 0, gate); }
void timeit_quest_5_0() { quest_two_qubit_gate(sv, 5, 0, gate); }
void timeit_quest_6_0() { quest_two_qubit_gate(sv, 6, 0, gate); }
void timeit_quest_7_0() { quest_two_qubit_gate(sv, 7, 0, gate); }
void timeit_quest_8_0() { quest_two_qubit_gate(sv, 8, 0, gate); }
void timeit_quest_9_0() { quest_two_qubit_gate(sv, 9, 0, gate); }
void timeit_quest_2_1() { quest_two_qubit_gate(sv, 2, 1, gate); }
void timeit_quest_3_1() { quest_two_qubit_gate(sv, 3, 1, gate); }
void timeit_quest_4_1() { quest_two_qubit_gate(sv, 4, 1, gate); }
void timeit_quest_5_1() { quest_two_qubit_gate(sv, 5, 1, gate); }
void timeit_quest_6_1() { quest_two_qubit_gate(sv, 6, 1, gate); }
void timeit_quest_7_1() { quest_two_qubit_gate(sv, 7, 1, gate); }
void timeit_quest_8_1() { quest_two_qubit_gate(sv, 8, 1, gate); }
void timeit_quest_9_1() { quest_two_qubit_gate(sv, 9, 1, gate); }
void timeit_quest_3_2() { quest_two_qubit_gate(sv, 3, 2, gate); }
void timeit_quest_4_2() { quest_two_qubit_gate(sv, 4, 2, gate); }
void timeit_quest_5_2() { quest_two_qubit_gate(sv, 5, 2, gate); }
void timeit_quest_6_2() { quest_two_qubit_gate(sv, 6, 2, gate); }
void timeit_quest_7_2() { quest_two_qubit_gate(sv, 7, 2, gate); }
void timeit_quest_8_2() { quest_two_qubit_gate(sv, 8, 2, gate); }
void timeit_quest_9_2() { quest_two_qubit_gate(sv, 9, 2, gate); }
void timeit_quest_4_3() { quest_two_qubit_gate(sv, 4, 3, gate); }
void timeit_quest_5_3() { quest_two_qubit_gate(sv, 5, 3, gate); }
void timeit_quest_6_3() { quest_two_qubit_gate(sv, 6, 3, gate); }
void timeit_quest_7_3() { quest_two_qubit_gate(sv, 7, 3, gate); }
void timeit_quest_8_3() { quest_two_qubit_gate(sv, 8, 3, gate); }
void timeit_quest_9_3() { quest_two_qubit_gate(sv, 9, 3, gate); }
void timeit_quest_5_4() { quest_two_qubit_gate(sv, 5, 4, gate); }
void timeit_quest_6_4() { quest_two_qubit_gate(sv, 6, 4, gate); }
void timeit_quest_7_4() { quest_two_qubit_gate(sv, 7, 4, gate); }
void timeit_quest_8_4() { quest_two_qubit_gate(sv, 8, 4, gate); }
void timeit_quest_9_4() { quest_two_qubit_gate(sv, 9, 4, gate); }
void timeit_quest_6_5() { quest_two_qubit_gate(sv, 6, 5, gate); }
void timeit_quest_7_5() { quest_two_qubit_gate(sv, 7, 5, gate); }
void timeit_quest_8_5() { quest_two_qubit_gate(sv, 8, 5, gate); }
void timeit_quest_9_5() { quest_two_qubit_gate(sv, 9, 5, gate); }
void timeit_quest_7_6() { quest_two_qubit_gate(sv, 7, 6, gate); }
void timeit_quest_8_6() { quest_two_qubit_gate(sv, 8, 6, gate); }
void timeit_quest_9_6() { quest_two_qubit_gate(sv, 9, 6, gate); }
void timeit_quest_8_7() { quest_two_qubit_gate(sv, 8, 7, gate); }
void timeit_quest_9_7() { quest_two_qubit_gate(sv, 9, 7, gate); }
void timeit_quest_9_8() { quest_two_qubit_gate(sv, 9, 8, gate); }


void (*timeit_quest[45])() = {
	timeit_quest_1_0,
	timeit_quest_2_0,
	timeit_quest_3_0,
	timeit_quest_4_0,
	timeit_quest_5_0,
	timeit_quest_6_0,
	timeit_quest_7_0,
	timeit_quest_8_0,
	timeit_quest_9_0,
	timeit_quest_2_1,
	timeit_quest_3_1,
	timeit_quest_4_1,
	timeit_quest_5_1,
	timeit_quest_6_1,
	timeit_quest_7_1,
	timeit_quest_8_1,
	timeit_quest_9_1,
	timeit_quest_3_2,
	timeit_quest_4_2,
	timeit_quest_5_2,
	timeit_quest_6_2,
	timeit_quest_7_2,
	timeit_quest_8_2,
	timeit_quest_9_2,
	timeit_quest_4_3,
	timeit_quest_5_3,
	timeit_quest_6_3,
	timeit_quest_7_3,
	timeit_quest_8_3,
	timeit_quest_9_3,
	timeit_quest_5_4,
	timeit_quest_6_4,
	timeit_quest_7_4,
	timeit_quest_8_4,
	timeit_quest_9_4,
	timeit_quest_6_5,
	timeit_quest_7_5,
	timeit_quest_8_5,
	timeit_quest_9_5,
	timeit_quest_7_6,
	timeit_quest_8_6,
	timeit_quest_9_6,
	timeit_quest_8_7,
	timeit_quest_9_7,
	timeit_quest_9_8,
};


int main(int argc, char *argv[]) {

	int num_threads = 1;
    if (argc > 1) {
        num_threads = atoi(argv[1]);
    }

	omp_set_num_threads(num_threads);

    TimingResult rst;

    printf("General:\n");
    for (int k = 0; k < 45; k++)
    {
		sleep(1);
        timeit(pre, timeit_general[k], post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    }

    printf("\n\nExplicit:\n");
    for (int k = 0; k < 45; k++)
    {
		sleep(1);
        timeit(pre, timeit_explicit[k], post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    }

    printf("\n\nQuest:\n");
    for (int k = 0; k < 45; k++)
    {
		sleep(1);
        timeit(pre, timeit_quest[k], post, &rst);
        printf("%f, ", rst.wall_stat[2] / rst.repetition * 1e6);
    }

    printf("\n");

}

