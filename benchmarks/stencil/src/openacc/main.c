
/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
#include <parboil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "file.h"
#include "common.h"
#include "kernels.h"

static int read_data ( float *A0, int nx, int ny, int nz, FILE *fp ) {
	int s = 0;
	for ( int i = 0; i < nz; ++i ) {
		for ( int j = 0; j < ny; ++j ) {
			for ( int k = 0; k < nx; ++k ) {
				fread ( A0 + s, sizeof ( float ), 1, fp );
				s++;
			}
		}
	}

	return 0;
}

int main ( int argc, char **argv ) {
	struct pb_TimerSet timers;
	struct pb_Parameters *parameters;

	printf ( "OPENACC-based 7 points stencil codes****\n" );
	parameters = pb_ReadParameters ( &argc, argv );

	pb_InitializeTimerSet ( &timers );
	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	//declaration
	float c0 = 1.0f / 6.0f;
	float c1 = 1.0f / 6.0f / 6.0f;

	if ( argc < 5 ) {
		printf ( 
			"Usage: probe nx ny nz tx ty t\n"
			"nx: the grid size x\n"
			"ny: the grid size y\n"
			"nz: the grid size z\n"
			"t: the iteration time\n"
		);

		return -1;
	}

	int nx = atoi ( argv[1] );
	if ( nx < 1 ) {
		return -1;
	}

	int ny = atoi ( argv[2] );
	if ( ny < 1 ) {
		return -1;
	}

	int nz = atoi ( argv[3] );
	if ( nz < 1 ) {
		return -1;
	}

	int iteration = atoi ( argv[4] );
	if ( iteration < 1 ) {
		return -1;
	}


	int size = nx * ny * nz;

	//host data
	float *h_A0 = ( float* )malloc ( sizeof ( float ) * size );
	float *h_Anext = ( float* )malloc ( sizeof ( float ) * size );
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );
	FILE *fp = fopen ( parameters->inpFiles[0], "rb" );
	read_data ( h_A0, nx, ny, nz, fp );
	fclose ( fp );

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
	memcpy ( h_Anext, h_A0, sizeof ( float ) * size );

	const unsigned int length = nx * ny * nz;
	#pragma acc data copyin(h_A0[0:length]) copy(h_Anext[0:length]) 
	{
		for ( int t = 0; t < iteration; ++t ) {
			cpu_stencil ( c0, c1, h_A0, h_Anext, nx, ny, nz );
			float *temp = h_A0;
			h_A0 = h_Anext;
			h_Anext = temp;
		}
	}

	if ( parameters->outFile ) {
		pb_SwitchToTimer ( &timers, pb_TimerID_IO );
		outputData ( parameters->outFile, h_Anext, nx, ny, nz );
	}
	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	free ( h_A0 );
	free ( h_Anext );
	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );

	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );

	return 0;
}
