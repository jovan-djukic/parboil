/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "computeQ.h"

void ComputePhiMagCPU ( int numK, float *phiR, float *phiI, float *phiMag ) {
    #pragma acc parallel loop present(phiR[0:numK], phiI[0:numK], phiMag[0:numK])
    for ( int indexK = 0; indexK < numK; ++indexK ) {
        float real = phiR[indexK];
        float imag = phiI[indexK];
        phiMag[indexK] = real * real + imag * imag;
    }
}

void ComputeQCPU ( int numK, int numX, struct kValues *kVals, float *x, float *y, float *z, float *Qr, float *Qi ) {
    #pragma acc parallel loop present(x[0:numX], y[0:numX], z[0:numX], Qr[0:numX], Qi[0:numX], kVals[0:numK]) 
    for ( int indexX = 0; indexX < numX; ++indexX ) {
        for ( int indexK = 0; indexK < numK; ++indexK ) {
            float expArg = PIx2 * ( kVals[indexK].Kx * x[indexX] + kVals[indexK].Ky * y[indexX] + kVals[indexK].Kz * z[indexX] );

            float cosArg = cosf ( expArg );
            float sinArg = sinf ( expArg );

            float phi = kVals[indexK].PhiMag;

            // #pragma acc atomic update
            Qr[indexX] += phi * cosArg;

            // #pragma acc atomic update
            Qi[indexX] += phi * sinArg;
        }
    }
}

void createDataStructsCPU ( int numK, int numX, float **phiMag, float **Qr, float **Qi ) {
    *phiMag = (float*) memalign ( 16, numK * sizeof ( float ) );
    *Qr = (float*) memalign ( 16, numX * sizeof ( float ) );
    *Qi = (float*) memalign ( 16, numX * sizeof ( float ) );

    // memset ( (void*) *Qr, 0, numX * sizeof ( float ) );
    // memset ( (void*) *Qi, 0, numX * sizeof ( float ) );
}
