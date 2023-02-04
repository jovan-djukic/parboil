/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/*############################################################################*/

#include "main.h"
#include "lbm.h"
#include <stdio.h>
#include <stdlib.h>

#include <sys/stat.h>

/*############################################################################*/
static LBM_Grid srcGrid, dstGrid;


/*############################################################################*/

#define REAL_MARGIN (CALC_INDEX(0, 0, 2, 0) - CALC_INDEX(0,0,0,0))
#define TOTAL_MARGIN (2*PADDED_X*PADDED_Y*N_CELL_ENTRIES)
#define TOTAL_SIZE (TOTAL_PADDED_CELLS*N_CELL_ENTRIES + 2*TOTAL_MARGIN)

struct pb_TimerSet timers;
int main( int nArgs, char* arg[] ) {
	MAIN_Param param;
	int t;

	pb_InitializeTimerSet(&timers);
    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);
    struct pb_Parameters* params;
    params = pb_ReadParameters(&nArgs, arg);
        

	MAIN_parseCommandLine( nArgs, arg, &param, params );
	MAIN_printInfo( &param );
	MAIN_initialize( &param );

	#pragma acc data copy(srcGrid[-REAL_MARGIN:TOTAL_SIZE]) copyin(dstGrid[-REAL_MARGIN:TOTAL_SIZE])
	{
		for( t = 1; t <= param.nTimeSteps; t++ ) {
			LBM_performStreamCollide( srcGrid, dstGrid );
			LBM_swapGrids( &srcGrid, &dstGrid );

			if( (t & 63) == 0 ) {
				printf( "timestep: %i\n", t );
	#if 0
				CUDA_LBM_getDeviceGrid((float**)&CUDA_srcGrid, (float**)&TEMP_srcGrid);
				LBM_showGridStatistics( *TEMP_srcGrid );
	#endif
			}
		}
	}

	MAIN_finalize( &param );

    pb_SwitchToTimer(&timers, pb_TimerID_NONE);
    pb_PrintTimerSet(&timers);
    pb_FreeParameters(params);
	return 0;
}

/*############################################################################*/

void MAIN_parseCommandLine( int nArgs, char* arg[], MAIN_Param* param, struct pb_Parameters * params ) {
	struct stat fileStat;

	if( nArgs < 2 ) {
		printf( "syntax: lbm <time steps>\n" );
		exit( 1 );
	}

	param->nTimeSteps     = atoi( arg[1] );

	if( params->inpFiles[0] != NULL ) {
		param->obstacleFilename = params->inpFiles[0];

		if( stat( param->obstacleFilename, &fileStat ) != 0 ) {
			printf( "MAIN_parseCommandLine: cannot stat obstacle file '%s'\n",
					param->obstacleFilename );
			exit( 1 );
		}
		if( fileStat.st_size != SIZE_X*SIZE_Y*SIZE_Z+(SIZE_Y+1)*SIZE_Z ) {
			printf( "MAIN_parseCommandLine:\n"
					"\tsize of file '%s' is %i bytes\n"
					"\texpected size is %i bytes\n",
					param->obstacleFilename, (int) fileStat.st_size,
					SIZE_X*SIZE_Y*SIZE_Z+(SIZE_Y+1)*SIZE_Z );
			exit( 1 );
		}
	}
	else param->obstacleFilename = NULL;

        param->resultFilename = params->outFile;
}

/*############################################################################*/

void MAIN_printInfo( const MAIN_Param* param ) {
	printf( "MAIN_printInfo:\n"
			"\tgrid size      : %i x %i x %i = %.2f * 10^6 Cells\n"
			"\tnTimeSteps     : %i\n"
			"\tresult file    : %s\n"
			"\taction         : %s\n"
			"\tsimulation type: %s\n"
			"\tobstacle file  : %s\n\n",
			SIZE_X, SIZE_Y, SIZE_Z, 1e-6*SIZE_X*SIZE_Y*SIZE_Z,
			param->nTimeSteps, param->resultFilename, 
			"store", "lid-driven cavity",
			(param->obstacleFilename == NULL) ? "<none>" :
			param->obstacleFilename );
}

/*############################################################################*/

void MAIN_initialize( const MAIN_Param* param ) {
	//Setup TEMP datastructures
	LBM_allocateGrid( (float**) &srcGrid );
	LBM_allocateGrid( (float**) &dstGrid );

	LBM_initializeGrid( srcGrid );
	LBM_initializeGrid( dstGrid );

	if( param->obstacleFilename != NULL ) {
    	pb_SwitchToTimer(&timers, pb_TimerID_IO);
		LBM_loadObstacleFile( srcGrid, param->obstacleFilename );
		LBM_loadObstacleFile( dstGrid, param->obstacleFilename );
    	pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);
	}

	LBM_initializeSpecialCellsForLDC( srcGrid );
	LBM_initializeSpecialCellsForLDC( dstGrid );

	LBM_showGridStatistics( srcGrid );
}

/*############################################################################*/

void MAIN_finalize( const MAIN_Param* param ) {
	LBM_showGridStatistics( srcGrid );

    pb_SwitchToTimer(&timers, pb_TimerID_IO);
	LBM_storeVelocityField( srcGrid, param->resultFilename, TRUE );
    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

	LBM_freeGrid( (float**) &srcGrid );
}

