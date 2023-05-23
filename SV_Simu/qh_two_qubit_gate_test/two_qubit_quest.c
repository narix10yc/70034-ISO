#include "two_qubit_quest.h"

void quest_two_qubit_gate(Statevector *sv, int k, int l, ComplexMatrix4 u) {

    // can't use qureg.stateVec as a private OMP var
    double *reVec = sv->real;
    double *imVec = sv->imag;
        
    size_t numTasks = sv->namp >> 2; // each iteration updates 4 amplitudes
    size_t thisTask;
    size_t thisGlobalInd00;
    size_t ind00, ind01, ind10, ind11;
    double re00, re01, re10, re11;
    double im00, im01, im10, im11;
    
    for (thisTask=0; thisTask<numTasks; thisTask++) {
        
        // determine ind00 of |..0..0..>
        ind00 = insertZeroBit(insertZeroBit(thisTask, l), k);
        ind01 = flipBit(ind00, k);
        ind10 = flipBit(ind00, l);
        ind11 = flipBit(ind01, l);

        // extract statevec amplitudes 
        re00 = reVec[ind00]; im00 = imVec[ind00];
        re01 = reVec[ind01]; im01 = imVec[ind01];
        re10 = reVec[ind10]; im10 = imVec[ind10];
        re11 = reVec[ind11]; im11 = imVec[ind11];

        // apply u * {amp00, amp01, amp10, amp11}
        reVec[ind00] = 
            u.real[0][0]*re00 - u.imag[0][0]*im00 +
            u.real[0][1]*re01 - u.imag[0][1]*im01 +
            u.real[0][2]*re10 - u.imag[0][2]*im10 +
            u.real[0][3]*re11 - u.imag[0][3]*im11;
        imVec[ind00] =
            u.imag[0][0]*re00 + u.real[0][0]*im00 +
            u.imag[0][1]*re01 + u.real[0][1]*im01 +
            u.imag[0][2]*re10 + u.real[0][2]*im10 +
            u.imag[0][3]*re11 + u.real[0][3]*im11;
            
        reVec[ind01] = 
            u.real[1][0]*re00 - u.imag[1][0]*im00 +
            u.real[1][1]*re01 - u.imag[1][1]*im01 +
            u.real[1][2]*re10 - u.imag[1][2]*im10 +
            u.real[1][3]*re11 - u.imag[1][3]*im11;
        imVec[ind01] =
            u.imag[1][0]*re00 + u.real[1][0]*im00 +
            u.imag[1][1]*re01 + u.real[1][1]*im01 +
            u.imag[1][2]*re10 + u.real[1][2]*im10 +
            u.imag[1][3]*re11 + u.real[1][3]*im11;
            
        reVec[ind10] = 
            u.real[2][0]*re00 - u.imag[2][0]*im00 +
            u.real[2][1]*re01 - u.imag[2][1]*im01 +
            u.real[2][2]*re10 - u.imag[2][2]*im10 +
            u.real[2][3]*re11 - u.imag[2][3]*im11;
        imVec[ind10] =
            u.imag[2][0]*re00 + u.real[2][0]*im00 +
            u.imag[2][1]*re01 + u.real[2][1]*im01 +
            u.imag[2][2]*re10 + u.real[2][2]*im10 +
            u.imag[2][3]*re11 + u.real[2][3]*im11;    
            
        reVec[ind11] = 
            u.real[3][0]*re00 - u.imag[3][0]*im00 +
            u.real[3][1]*re01 - u.imag[3][1]*im01 +
            u.real[3][2]*re10 - u.imag[3][2]*im10 +
            u.real[3][3]*re11 - u.imag[3][3]*im11;
        imVec[ind11] =
            u.imag[3][0]*re00 + u.real[3][0]*im00 +
            u.imag[3][1]*re01 + u.real[3][1]*im01 +
            u.imag[3][2]*re10 + u.real[3][2]*im10 +
            u.imag[3][3]*re11 + u.real[3][3]*im11;    
        
    }
}


