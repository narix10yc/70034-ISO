#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include "qh_type.h"

Statevector *svCreate(int nqubits);

void svDestroy(Statevector *sv);

void svInit(Statevector *sv);

void svRandomize(Statevector *sv);

void svPrintAmplitude(Statevector *sv);

void svPrintProbability(Statevector *sv);


#endif