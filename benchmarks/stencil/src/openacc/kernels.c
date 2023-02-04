/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include "common.h"

void cpu_stencil ( float c0, float c1, float *A0, float * Anext, const int nx, const int ny, const int nz ) {
	const unsigned int length = nx * ny * nz;
	const unsigned int iLimit = nx - 2;
	const unsigned int jLimit = ny - 2;
	const unsigned int kLimit = nz - 2;

	#pragma acc parallel loop collapse(3) present(A0[0:length], Anext[0:length])
	for ( unsigned int k = 0; k < kLimit; ++k ) {
		for ( unsigned int j = 0; j < jLimit; ++j ) {
			for ( unsigned int i = 0; i < iLimit; ++i ) {
				const unsigned int iIndex = i + 1;
				const unsigned int jIndex = j + 1;
				const unsigned int kIndex = k + 1;

				const float kUp    = A0[Index3D (nx, ny, iIndex    , jIndex    , kIndex + 1)];
				const float kDown  = A0[Index3D (nx, ny, iIndex    , jIndex    , kIndex - 1)];
				const float jUp    = A0[Index3D (nx, ny, iIndex    , jIndex + 1, kIndex    )];
				const float jDown  = A0[Index3D (nx, ny, iIndex    , jIndex - 1, kIndex    )];
				const float iUp    = A0[Index3D (nx, ny, iIndex + 1, jIndex    , kIndex    )];
				const float iDown  = A0[Index3D (nx, ny, iIndex - 1, jIndex    , kIndex    )];
				const float center = A0[Index3D (nx, ny, iIndex    , jIndex    , kIndex    )];

				Anext[Index3D (nx, ny, iIndex, jIndex, kIndex)] =  ( kUp + kDown + jUp + jDown + iUp + iDown ) * c1 - center * c0;
			}
		}
	}
}


