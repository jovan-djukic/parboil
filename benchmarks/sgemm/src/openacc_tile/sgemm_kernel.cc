/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/* 
 * Base C implementation of MM
 */

#include <stdlib.h>
#include <iostream>

void basicSgemm ( char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc ) {
    if ( ( transa != 'N' ) && ( transa != 'n' ) ) {
        std::cerr << "unsupported value of 'transa' in regtileSgemm ( )" << std::endl;
        return;
    }

    if ( ( transb != 'T' ) && ( transb != 't' ) ) {
        std::cerr << "unsupported value of 'transb' in regtileSgemm ( )" << std::endl;
        return;
    }

    const int aLength = m * k;
    const int bLength = k * n;
    const int cLength = m * n;

    #pragma acc data copyin(A [0:aLength], B [0:bLength]) copyout(C [0:cLength])
    {
        #pragma acc parallel loop tile(16, 16)
        for ( int mm = 0; mm < m; ++mm ) {
            for ( int nn = 0; nn < n; ++nn ) {
                float temp = 0;

                #pragma acc loop reduction(+: temp)
                for ( int i = 0; i < k; ++i ) {
                    temp += A[mm + i * lda] * B[nn + i * ldb];
                }

                C[mm + nn * ldc] = C[mm + nn * ldc] * beta + alpha * temp;
            }
        }
    }
}
