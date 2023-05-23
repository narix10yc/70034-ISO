#ifndef TWO_QUBIT_QUEST_H
#define TWO_QUBIT_QUEST_H

#include "qh_type.h"

static inline long long int flipBit(const long long int number, const int bitInd) {
    return (number ^ (1LL << bitInd));
}

static inline long long int insertZeroBit(const long long int number, const int index) {
    long long int left, right;
    left = (number >> index) << index;
    right = number - left;
    return (left << 1) ^ right;
}

void quest_two_qubit_gate(Statevector *sv, int k, int l, ComplexMatrix4 u);

#endif