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
    float *QrTemp = memalign ( 16, numK * numX * sizeof ( float ) );
    float *QiTemp = memalign ( 16, numK * numX * sizeof ( float ) );

    #pragma acc data present(x[0:numX], y[0:numX], z[0:numX], Qr[0:numX], Qi[0:numX], kVals[0:numK]) create(QrTemp[0:(numK * numX)], QiTemp[0:(numK * numX)])
    {
        #pragma acc parallel loop collapse(2)
        for ( int indexX = 0; indexX < numX; ++indexX ) {
            for ( int indexK = 0; indexK < numK; ++indexK ) {
                float expArg = PIx2 * ( kVals[indexK].Kx * x[indexX] + kVals[indexK].Ky * y[indexX] + kVals[indexK].Kz * z[indexX] );

                float cosArg = cosf ( expArg );
                float sinArg = sinf ( expArg );

                float phi = kVals[indexK].PhiMag;

                int index = indexK * numX + indexX;

                QrTemp[index] = phi * cosArg;
                QiTemp[index] = phi * sinArg;
            }
        }

        #pragma acc parallel loop 
        for ( int indexX = 0; indexX < numX; ++indexX ) {
            float QrTotal = 0;
            float QiTotal = 0;
            #pragma acc loop reduction(+:QrTotal,QiTotal)
            for ( int indexK = 0; indexK < numK; ++indexK ) {
                QrTotal += QrTemp[indexK * numX + indexX];
                QiTotal += QiTemp[indexK * numX + indexX];
            }

            Qr[indexX] = QrTotal;
            Qi[indexX] = QiTotal;
        }
    }

    free ( QrTemp );
    free ( QiTemp );
}

void createDataStructsCPU ( int numK, int numX, float **phiMag, float **Qr, float **Qi ) {
    *phiMag = (float*) memalign ( 16, numK * sizeof ( float ) );
    *Qr = (float*) memalign ( 16, numX * sizeof ( float ) );
    *Qi = (float*) memalign ( 16, numX * sizeof ( float ) );

    // memset ( (void*) *Qr, 0, numX * sizeof ( float ) );
    // memset ( (void*) *Qi, 0, numX * sizeof ( float ) );
}
