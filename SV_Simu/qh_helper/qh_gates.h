#ifndef QH_GATES_H
#define QH_GATES_H

#include "qh_type.h"

void applyGateSingleUnitary(Statevector *sv, int k, ComplexMatrix2 u);

void applyGateX(Statevector *sv, int k);

#endif
