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

#define ROW_WISE_INDEX(i, j, rows, columns) ( (i) * (columns) + (j) )
#define COLUMN_WISE_INDEX(i, j, rows, columns) ( (i) + (j) * (rows) )
#define TCOLUMN_WISE_INDEX(i, j, rows, columns) ROW_WISE_INDEX(i, j, rows, columns)

#define B_BLOCK_ROWS 8
#define C_BLOCK_LENGTH 16

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

    #pragma acc parallel loop collapse(2) copyin(A[0:aLength], B[0:bLength]) copyout(C[0:cLength])
    for ( int tx = 0; tx < m; ++tx ) {
        for ( int ty = 0; ty < n; ty += C_BLOCK_LENGTH ) {
            float cTemp[C_BLOCK_LENGTH] = { 0 };

            for ( int bx = 0; bx < k; bx += B_BLOCK_ROWS ) {
                const int sharedStart0 = COLUMN_WISE_INDEX(bx + 0, ty + 0, k, n);
                const int sharedStart1 = COLUMN_WISE_INDEX(bx + 1, ty + 0, k, n);
                const int sharedStart2 = COLUMN_WISE_INDEX(bx + 2, ty + 0, k, n);
                const int sharedStart3 = COLUMN_WISE_INDEX(bx + 3, ty + 0, k, n);
                const int sharedStart4 = COLUMN_WISE_INDEX(bx + 4, ty + 0, k, n);
                const int sharedStart5 = COLUMN_WISE_INDEX(bx + 5, ty + 0, k, n);
                const int sharedStart6 = COLUMN_WISE_INDEX(bx + 6, ty + 0, k, n);
                const int sharedStart7 = COLUMN_WISE_INDEX(bx + 7, ty + 0, k, n);

                #pragma acc cache(B [sharedStart0 : C_BLOCK_LENGTH]\
                                    [sharedStart1 : C_BLOCK_LENGTH]\
                                    [sharedStart2 : C_BLOCK_LENGTH]\
                                    [sharedStart3 : C_BLOCK_LENGTH]\
                                    [sharedStart4 : C_BLOCK_LENGTH]\
                                    [sharedStart5 : C_BLOCK_LENGTH]\
                                    [sharedStart6 : C_BLOCK_LENGTH]\
                                    [sharedStart7 : C_BLOCK_LENGTH]\
                )
                for ( int i = 0; i < B_BLOCK_ROWS; ++i ) { 
                    float a = A[COLUMN_WISE_INDEX(tx, bx + i, m, k)];
                    for ( int j = 0; j < C_BLOCK_LENGTH; ++j ) {
                        cTemp[j] += a * B[TCOLUMN_WISE_INDEX(bx + i, ty + j, k, n)];
                    }
                }

            }

            for ( int i = 0; i < C_BLOCK_LENGTH; ++i ) {
                C[COLUMN_WISE_INDEX(tx, ty + i, m, n)] = C[COLUMN_WISE_INDEX(tx, ty + i, m, n)] * beta + cTemp[i] * alpha;
            }
        }
    }
}
