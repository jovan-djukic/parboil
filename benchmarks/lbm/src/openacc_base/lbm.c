/* $Id: lbm.c,v 1.1 2008/03/04 17:30:02 stratton Exp $ */

/*############################################################################*/

#include "lbm.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

/*############################################################################*/

#define DFL1 ( 1.0 /  3.0 )
#define DFL2 ( 1.0 / 18.0 )
#define DFL3 ( 1.0 / 36.0 )

/*############################################################################*/

void LBM_allocateGrid ( float **ptr ) {
	const size_t margin = 2 * SIZE_X * SIZE_Y * N_CELL_ENTRIES;
	// const size_t size   = sizeof ( LBM_Grid ) + 2 * margin * sizeof ( float );
	const size_t size   = SIZE_X * SIZE_Y * SIZE_Z * N_CELL_ENTRIES * sizeof ( float ) + 2 * margin * sizeof ( float );

	*ptr = malloc ( size );
	if ( !(*ptr) ) {
		printf ("LBM_allocateGrid: could not allocate %.1f MByte\n", size / ( 1024.0 * 1024.0 ) );
		exit ( 1 );
	}
#if !defined ( SPEC_CPU )
	printf ( "LBM_allocateGrid: allocated %.1f MByte\n", size / ( 1024.0 * 1024.0 ) );
#endif

	*ptr += margin;
}

/*############################################################################*/

void LBM_freeGrid ( float **ptr ) {
	const size_t margin = 2 * SIZE_X * SIZE_Y * N_CELL_ENTRIES;

	free ( *ptr - margin );
	*ptr = NULL;
}

/*############################################################################*/

// void LBM_initializeGrid ( LBM_Grid grid ) {
void LBM_initializeGrid ( float *grid ) {
	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for
#endif
#endif
	SWEEP_START ( 0, 0, -2, 0, 0, SIZE_Z + 2 )
		LOCAL ( grid, C  ) = DFL1;
		LOCAL ( grid, N  ) = DFL2;
		LOCAL ( grid, S  ) = DFL2;
		LOCAL ( grid, E  ) = DFL2;
		LOCAL ( grid, W  ) = DFL2;
		LOCAL ( grid, T  ) = DFL2;
		LOCAL ( grid, B  ) = DFL2;
		LOCAL ( grid, NE ) = DFL3;
		LOCAL ( grid, NW ) = DFL3;
		LOCAL ( grid, SE ) = DFL3;
		LOCAL ( grid, SW ) = DFL3;
		LOCAL ( grid, NT ) = DFL3;
		LOCAL ( grid, NB ) = DFL3;
		LOCAL ( grid, ST ) = DFL3;
		LOCAL ( grid, SB ) = DFL3;
		LOCAL ( grid, ET ) = DFL3;
		LOCAL ( grid, EB ) = DFL3;
		LOCAL ( grid, WT ) = DFL3;
		LOCAL ( grid, WB ) = DFL3;

		CLEAR_ALL_FLAGS_SWEEP ( grid );
	SWEEP_END
}

/*############################################################################*/

// void LBM_swapGrids ( LBM_GridPtr *grid1, LBM_GridPtr *grid2 ) {
void LBM_swapGrids ( float **grid1, float **grid2 ) {
	float *auxillary = *grid1;
	*grid1 = *grid2;
	*grid2 = auxillary;
}

/*############################################################################*/

// void LBM_loadObstacleFile ( LBM_Grid grid, const char *filename ) {
void LBM_loadObstacleFile ( float *grid, const char *filename ) {
	FILE *file = fopen(filename, "rb");

	for ( int z = 0; z < SIZE_Z; ++z ) {
		for ( int y = 0; y < SIZE_Y; ++y ) {
			for ( int x = 0; x < SIZE_X; ++x ) {
				if ( fgetc ( file ) != '.' ) {
					SET_FLAG ( grid, x, y, z, OBSTACLE );
				}
			}
			fgetc ( file );
		}
		fgetc ( file );
	}

	fclose ( file );
}

/*############################################################################*/

// void LBM_initializeSpecialCellsForLDC ( LBM_Grid grid ) {
void LBM_initializeSpecialCellsForLDC ( float *grid ) {
	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for private(x, y)
#endif
#endif
	for ( int z = -2; z < ( SIZE_Z + 2 ); ++z ) {
		for ( int y = 0; y < SIZE_Y; ++y ) {
			for ( int x = 0; x < SIZE_X; ++x ) {
				if ( x == 0 || x == ( SIZE_X - 1 ) || y == 0 || y == ( SIZE_Y - 1 ) || z == 0 || z == ( SIZE_Z - 1 ) ) {
					SET_FLAG ( grid, x, y, z, OBSTACLE );
				} else {
					if ( ( z == 1 || z == ( SIZE_Z - 2 ) ) && x > 1 && x < ( SIZE_X - 2 ) && y > 1 && y < ( SIZE_Y - 2 ) ) {
						SET_FLAG ( grid, x, y, z, ACCEL );
					}
				}
			}
		}
	}
}

/*############################################################################*/

// void LBM_initializeSpecialCellsForChannel ( LBM_Grid grid ) {
void LBM_initializeSpecialCellsForChannel ( float *grid ) {
	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for private(x, y)
#endif
#endif
	for ( int z = -2; z < ( SIZE_Z + 2 ); ++z ) {
		for ( int y = 0; y < SIZE_Y; ++y ) {
			for ( int x = 0; x < SIZE_X; ++x ) {
				if ( x == 0 || x == ( SIZE_X - 1 ) || y == 0 || y == ( SIZE_Y - 1 ) ) {
					SET_FLAG ( grid, x, y, z, OBSTACLE );

					if ( ( z == 0 || z == ( SIZE_Z - 1 ) ) && !TEST_FLAG ( grid, x, y, z, OBSTACLE ) ) {
						SET_FLAG ( grid, x, y, z, IN_OUT_FLOW );
					}
				}
			}
		}
	}
}

/*############################################################################*/

// void LBM_performStreamCollide ( LBM_Grid srcGrid, LBM_Grid dstGrid ) {
void LBM_performStreamCollide ( float *srcGrid, float *dstGrid ) {

	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for private(ux, uy, uz, u2, rho)
#endif
#endif

	const size_t margin = 2 * SIZE_X * SIZE_Y * N_CELL_ENTRIES;
	const size_t size   = SIZE_X * SIZE_Y * SIZE_Z * N_CELL_ENTRIES;

	// TODO: play with vectors and gangs
	#pragma acc parallel loop present(srcGrid[-margin:(size + 2 * margin)], dstGrid[-margin:(size + 2 * margin)])
	SWEEP_START(0, 0, 0, 0, 0, SIZE_Z)
		float srcC  = SRC_C  ( srcGrid );
		float srcN  = SRC_N  ( srcGrid );
		float srcS  = SRC_S  ( srcGrid );
		float srcE  = SRC_E  ( srcGrid );
		float srcW  = SRC_W  ( srcGrid );
		float srcT  = SRC_T  ( srcGrid );
		float srcB  = SRC_B  ( srcGrid );
		float srcNE = SRC_NE ( srcGrid );
		float srcNW = SRC_NW ( srcGrid );
		float srcSE = SRC_SE ( srcGrid );
		float srcSW = SRC_SW ( srcGrid );
		float srcNT = SRC_NT ( srcGrid );
		float srcNB = SRC_NB ( srcGrid );
		float srcST = SRC_ST ( srcGrid );
		float srcSB = SRC_SB ( srcGrid );
		float srcET = SRC_ET ( srcGrid );
		float srcEB = SRC_EB ( srcGrid );
		float srcWT = SRC_WT ( srcGrid );
		float srcWB = SRC_WB ( srcGrid );

		//Test whether the cell is fluid or obstacle
		if ( TEST_FLAG_SWEEP ( srcGrid, OBSTACLE ) ) {
			//Swizzle the inputs: reflect any fluid coming into this cell back to where it came from
			float swap;
			swap = srcN ; srcN  = srcS ; srcS  = swap;
			swap = srcE ; srcE  = srcW ; srcW  = swap;
			swap = srcT ; srcT  = srcB ; srcB  = swap;
			swap = srcNE; srcNE = srcSW; srcSW = swap;
			swap = srcNW; srcNW = srcSE; srcSE = swap;
			swap = srcNT; srcNT = srcSB; srcSB = swap; 
			swap = srcNB; srcNB = srcST; srcST = swap;
			swap = srcET; srcET = srcWB; srcWB = swap;
			swap = srcEB; srcEB = srcWT; srcWT = swap;
		} else {
			float rho =   srcC  + srcN
				        + srcS  + srcE
				        + srcW  + srcT
				        + srcB  + srcNE
				        + srcNW + srcSE
				        + srcSW + srcNT
				        + srcNB + srcST
				        + srcSB + srcET
				        + srcEB + srcWT
				        + srcWB;

			float ux = 	+ srcE  - srcW
						+ srcNE - srcNW
						+ srcSE - srcSW
						+ srcET + srcEB
						- srcWT - srcWB;

			float uy =  + srcN  - srcS
						+ srcNE + srcNW
						- srcSE - srcSW
						+ srcNT + srcNB
						- srcST - srcSB;

			float uz =  + srcT  - srcB
						+ srcNT - srcNB
						+ srcST - srcSB
						+ srcET - srcEB
						+ srcWT - srcWB;

			ux /= rho;
			uy /= rho;
			uz /= rho;

			if ( TEST_FLAG_SWEEP ( srcGrid, ACCEL ) ) {
				ux = 0.005f;
				uy = 0.002f;
				uz = 0.000f;
			}

			float u2 = 1.5f * ( ux * ux + uy * uy + uz * uz ) - 1.0f;

			float base = OMEGA * rho;
			float dfl1 = DFL1 * base;
		    float dfl2 = DFL2 * base;	
			float dfl3 = DFL3 * base;

			float oneMinusOmega = 1.0f-OMEGA;


			//Put the output values for this cell in the shared memory
			srcC  = oneMinusOmega * srcC  + dfl1 * (                                                 - u2 );
			srcN  = oneMinusOmega * srcN  + dfl2 * ( (        uy ) * ( 4.5f * (        uy ) + 3.0f ) - u2 );
			srcS  = oneMinusOmega * srcS  + dfl2 * ( (        uy ) * ( 4.5f * (        uy ) - 3.0f ) - u2 );
			srcT  = oneMinusOmega * srcT  + dfl2 * ( (        uz ) * ( 4.5f * (        uz ) + 3.0f ) - u2 );
			srcB  = oneMinusOmega * srcB  + dfl2 * ( (        uz ) * ( 4.5f * (        uz ) - 3.0f ) - u2 );
			srcE  = oneMinusOmega * srcE  + dfl2 * ( (        ux ) * ( 4.5f * (        ux ) + 3.0f ) - u2 );
			srcW  = oneMinusOmega * srcW  + dfl2 * ( (        ux ) * ( 4.5f * (        ux ) - 3.0f ) - u2 );
			srcNT = oneMinusOmega * srcNT + dfl3 * ( ( + uy + uz ) * ( 4.5f * ( + uy + uz ) + 3.0f ) - u2 );
			srcNB = oneMinusOmega * srcNB + dfl3 * ( ( + uy - uz ) * ( 4.5f * ( + uy - uz ) + 3.0f ) - u2 );
			srcST = oneMinusOmega * srcST + dfl3 * ( ( - uy + uz ) * ( 4.5f * ( - uy + uz ) + 3.0f ) - u2 );
			srcSB = oneMinusOmega * srcSB + dfl3 * ( ( - uy - uz ) * ( 4.5f * ( - uy - uz ) + 3.0f ) - u2 );
			srcNE = oneMinusOmega * srcNE + dfl3 * ( ( + ux + uy ) * ( 4.5f * ( + ux + uy ) + 3.0f ) - u2 );
			srcSE = oneMinusOmega * srcSE + dfl3 * ( ( + ux - uy ) * ( 4.5f * ( + ux - uy ) + 3.0f ) - u2 );
			srcET = oneMinusOmega * srcET + dfl3 * ( ( + ux + uz ) * ( 4.5f * ( + ux + uz ) + 3.0f ) - u2 );
			srcEB = oneMinusOmega * srcEB + dfl3 * ( ( + ux - uz ) * ( 4.5f * ( + ux - uz ) + 3.0f ) - u2 );
			srcNW = oneMinusOmega * srcNW + dfl3 * ( ( - ux + uy ) * ( 4.5f * ( - ux + uy ) + 3.0f ) - u2 );
			srcSW = oneMinusOmega * srcSW + dfl3 * ( ( - ux - uy ) * ( 4.5f * ( - ux - uy ) + 3.0f ) - u2 );
			srcWT = oneMinusOmega * srcWT + dfl3 * ( ( - ux + uz ) * ( 4.5f * ( - ux + uz ) + 3.0f ) - u2 );
			srcWB = oneMinusOmega * srcWB + dfl3 * ( ( - ux - uz ) * ( 4.5f * ( - ux - uz ) + 3.0f ) - u2 );
		}

		// Write the results computed above
		//This is a scatter operation of the SCATTER preprocessor variable is defined in layout_config.h, or a "local" write otherwise
		DST_C  ( dstGrid ) = srcC;
		DST_N  ( dstGrid ) = srcN; 
		DST_S  ( dstGrid ) = srcS;
		DST_E  ( dstGrid ) = srcE;
		DST_W  ( dstGrid ) = srcW;
		DST_T  ( dstGrid ) = srcT;
		DST_B  ( dstGrid ) = srcB;
		DST_NE ( dstGrid ) = srcNE;
		DST_NW ( dstGrid ) = srcNW;
		DST_SE ( dstGrid ) = srcSE;
		DST_SW ( dstGrid ) = srcSW;
		DST_NT ( dstGrid ) = srcNT;
		DST_NB ( dstGrid ) = srcNB;
		DST_ST ( dstGrid ) = srcST;
		DST_SB ( dstGrid ) = srcSB;
		DST_ET ( dstGrid ) = srcET;
		DST_EB ( dstGrid ) = srcEB;
		DST_WT ( dstGrid ) = srcWT;
		DST_WB ( dstGrid ) = srcWB;
	SWEEP_END
}

/*############################################################################*/

// void LBM_handleInOutFlow ( LBM_Grid srcGrid ) {
void LBM_handleInOutFlow ( float *srcGrid ) {
	/* inflow */
	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for private(ux, uy, uz, rho, ux1, uy1, uz1, rho1, ux2, uy2, uz2, rho2, u2, px, py)
#endif
#endif
	SWEEP_START ( 0, 0, 0, 0, 0, 1 ) 
		float rho1 = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, C ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, N ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, E ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, T ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, SE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, ST ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, WT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 1, WB );
		float rho2 = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, C ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, N ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, E ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, T ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, SE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, ST ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, WT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, 2, WB );

		float rho = 2.0 * rho1 - rho2;

		float px = ( SWEEP_X / ( 0.5 * ( SIZE_X - 1 ) ) ) - 1.0;
		float py = ( SWEEP_Y / ( 0.5 * ( SIZE_Y - 1 ) ) ) - 1.0;
		float ux = 0.00;
		float uy = 0.00;
		float uz = 0.01 * ( 1.0 - px * px ) * ( 1.0 - py * py );

		float u2 = 1.5 * ( ux * ux + uy * uy + uz * uz );

		LOCAL ( srcGrid, C  ) = DFL1 * rho * ( 1.0 - u2 );

		LOCAL ( srcGrid, N  ) = DFL2 * rho * ( 1.0 + uy * ( 4.5 * uy + 3.0 ) - u2 );
		LOCAL ( srcGrid, S  ) = DFL2 * rho * ( 1.0 + uy * ( 4.5 * uy - 3.0 ) - u2 );
		LOCAL ( srcGrid, E  ) = DFL2 * rho * ( 1.0 + ux * ( 4.5 * ux + 3.0 ) - u2 );
		LOCAL ( srcGrid, W  ) = DFL2 * rho * ( 1.0 + ux * ( 4.5 * ux - 3.0 ) - u2 );
		LOCAL ( srcGrid, T  ) = DFL2 * rho * ( 1.0 + uz * ( 4.5 * uz + 3.0 ) - u2 );
		LOCAL ( srcGrid, B  ) = DFL2 * rho * ( 1.0 + uz * ( 4.5 * uz - 3.0 ) - u2 );

		LOCAL ( srcGrid, NE ) = DFL3 * rho * ( 1.0 + ( +ux + uy ) * ( 4.5 * ( +ux + uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NW ) = DFL3 * rho * ( 1.0 + ( -ux + uy ) * ( 4.5 * ( -ux + uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SE ) = DFL3 * rho * ( 1.0 + ( +ux - uy ) * ( 4.5 * ( +ux - uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SW ) = DFL3 * rho * ( 1.0 + ( -ux - uy ) * ( 4.5 * ( -ux - uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NT ) = DFL3 * rho * ( 1.0 + ( +uy + uz ) * ( 4.5 * ( +uy + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NB ) = DFL3 * rho * ( 1.0 + ( +uy - uz ) * ( 4.5 * ( +uy - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, ST ) = DFL3 * rho * ( 1.0 + ( -uy + uz ) * ( 4.5 * ( -uy + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SB ) = DFL3 * rho * ( 1.0 + ( -uy - uz ) * ( 4.5 * ( -uy - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, ET ) = DFL3 * rho * ( 1.0 + ( +ux + uz ) * ( 4.5 * ( +ux + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, EB ) = DFL3 * rho * ( 1.0 + ( +ux - uz ) * ( 4.5 * ( +ux - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, WT ) = DFL3 * rho * ( 1.0 + ( -ux + uz ) * ( 4.5 * ( -ux + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, WB ) = DFL3 * rho * ( 1.0 + ( -ux - uz ) * ( 4.5 * ( -ux - uz ) + 3.0 ) - u2 );
	SWEEP_END

	/* outflow */
	/*voption indep*/
#if !defined ( SPEC_CPU )
#ifdef _OPENMP
#pragma omp parallel for private(ux, uy, uz, rho, ux1, uy1, uz1, rho1, ux2, uy2, uz2, rho2, u2, px, py)
#endif
#endif

	SWEEP_START ( 0, 0, SIZE_Z - 1, 0, 0, SIZE_Z)
		float rho1 = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, C ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, N ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, E ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, T ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ST ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WB );
		float ux1  = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, E ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, EB ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WB );
		float uy1  = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, N ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NW ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NB ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ST ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SB );
		float uz1  = GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, T ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ST ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, ET ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -1, WB );

		ux1 /= rho1;
		uy1 /= rho1;
		uz1 /= rho1;

		float rho2 = +GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, C ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, N ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, E ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, T ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ST ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WB );
		float ux2  = +GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, E ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, W ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ET ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, EB ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WB );
		float uy2  = +GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, N ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, S ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NE ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NW ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SE ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SW ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NT ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NB ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ST ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SB );
		float uz2  = +GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, T ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, B ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, NB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ST ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, SB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, ET ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, EB ) + GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WT ) - GRID_ENTRY_SWEEP ( srcGrid, 0, 0, -2, WB );

		ux2 /= rho2;
		uy2 /= rho2;
		uz2 /= rho2;

		float rho = 1.0;

		float ux = 2 * ux1 - ux2;
		float uy = 2 * uy1 - uy2;
		float uz = 2 * uz1 - uz2;

		float u2 = 1.5 * ( ux * ux + uy * uy + uz * uz );

		LOCAL ( srcGrid, C ) = DFL1 * rho * ( 1.0 - u2 );

		LOCAL ( srcGrid, N ) = DFL2 * rho * ( 1.0 + uy * ( 4.5 * uy + 3.0 ) - u2 );
		LOCAL ( srcGrid, S ) = DFL2 * rho * ( 1.0 + uy * ( 4.5 * uy - 3.0 ) - u2 );
		LOCAL ( srcGrid, E ) = DFL2 * rho * ( 1.0 + ux * ( 4.5 * ux + 3.0 ) - u2 );
		LOCAL ( srcGrid, W ) = DFL2 * rho * ( 1.0 + ux * ( 4.5 * ux - 3.0 ) - u2 );
		LOCAL ( srcGrid, T ) = DFL2 * rho * ( 1.0 + uz * ( 4.5 * uz + 3.0 ) - u2 );
		LOCAL ( srcGrid, B ) = DFL2 * rho * ( 1.0 + uz * ( 4.5 * uz - 3.0 ) - u2 );

		LOCAL ( srcGrid, NE ) = DFL3 * rho * ( 1.0 + ( +ux + uy ) * ( 4.5 * ( +ux + uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NW ) = DFL3 * rho * ( 1.0 + ( -ux + uy ) * ( 4.5 * ( -ux + uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SE ) = DFL3 * rho * ( 1.0 + ( +ux - uy ) * ( 4.5 * ( +ux - uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SW ) = DFL3 * rho * ( 1.0 + ( -ux - uy ) * ( 4.5 * ( -ux - uy ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NT ) = DFL3 * rho * ( 1.0 + ( +uy + uz ) * ( 4.5 * ( +uy + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, NB ) = DFL3 * rho * ( 1.0 + ( +uy - uz ) * ( 4.5 * ( +uy - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, ST ) = DFL3 * rho * ( 1.0 + ( -uy + uz ) * ( 4.5 * ( -uy + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, SB ) = DFL3 * rho * ( 1.0 + ( -uy - uz ) * ( 4.5 * ( -uy - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, ET ) = DFL3 * rho * ( 1.0 + ( +ux + uz ) * ( 4.5 * ( +ux + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, EB ) = DFL3 * rho * ( 1.0 + ( +ux - uz ) * ( 4.5 * ( +ux - uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, WT ) = DFL3 * rho * ( 1.0 + ( -ux + uz ) * ( 4.5 * ( -ux + uz ) + 3.0 ) - u2 );
		LOCAL ( srcGrid, WB ) = DFL3 * rho * ( 1.0 + ( -ux - uz ) * ( 4.5 * ( -ux - uz ) + 3.0 ) - u2 );
	SWEEP_END
}

/*############################################################################*/

// void LBM_showGridStatistics ( LBM_Grid grid ) {
void LBM_showGridStatistics ( float *grid ) {
	int nObstacleCells = 0;
	int	nAccelCells = 0;
	int	nFluidCells = 0;
	float minU2 = 1e+30, maxU2 = -1e+30;
	float minRho = 1e+30, maxRho = -1e+30;
	float mass = 0;

	SWEEP_START ( 0, 0, 0, 0, 0, SIZE_Z ) 
		float rho = LOCAL ( grid, C ) + LOCAL ( grid, N ) + LOCAL ( grid, S ) + LOCAL ( grid, E ) + LOCAL ( grid, W ) + LOCAL ( grid, T ) + LOCAL ( grid, B ) + LOCAL ( grid, NE ) + LOCAL ( grid, NW ) + LOCAL ( grid, SE ) + LOCAL ( grid, SW ) + LOCAL ( grid, NT ) + LOCAL ( grid, NB ) + LOCAL ( grid, ST ) + LOCAL ( grid, SB ) + LOCAL ( grid, ET ) + LOCAL ( grid, EB ) + LOCAL ( grid, WT ) + LOCAL ( grid, WB );

		if ( rho < minRho ) {
			minRho = rho;
		}
		if ( rho > maxRho ) {
			maxRho = rho;
		}

		mass += rho;

		if ( TEST_FLAG_SWEEP ( grid, OBSTACLE ) ) {
			nObstacleCells++;
		} else {
			if ( TEST_FLAG_SWEEP ( grid, ACCEL ) ) {
				nAccelCells++;
			} else {
				nFluidCells++;
			}

			float ux = LOCAL ( grid, E ) - LOCAL ( grid, W ) + LOCAL ( grid, NE ) - LOCAL ( grid, NW ) + LOCAL ( grid, SE ) - LOCAL ( grid, SW ) + LOCAL ( grid, ET ) + LOCAL ( grid, EB ) - LOCAL ( grid, WT ) - LOCAL ( grid, WB );
			float uy = LOCAL ( grid, N ) - LOCAL ( grid, S ) + LOCAL ( grid, NE ) + LOCAL ( grid, NW ) - LOCAL ( grid, SE ) - LOCAL ( grid, SW ) + LOCAL ( grid, NT ) + LOCAL ( grid, NB ) - LOCAL ( grid, ST ) - LOCAL ( grid, SB );
			float uz = LOCAL ( grid, T ) - LOCAL ( grid, B ) + LOCAL ( grid, NT ) - LOCAL ( grid, NB ) + LOCAL ( grid, ST ) - LOCAL ( grid, SB ) + LOCAL ( grid, ET ) - LOCAL ( grid, EB ) + LOCAL ( grid, WT ) - LOCAL ( grid, WB );
			float u2 = ( ux * ux + uy * uy + uz * uz ) / ( rho * rho );

			if ( u2 < minU2 ) {
				minU2 = u2;
			}
			if ( u2 > maxU2 ) {
				maxU2 = u2;
			}
		}
	SWEEP_END

	printf ( 
		"LBM_showGridStatistics:\n"
			"\tnObstacleCells: %7i nAccelCells: %7i nFluidCells: %7i\n"
			"\tminRho: %8.4f maxRho: %8.4f mass: %e\n"
			"\tminU: %e maxU: %e\n\n",
		   nObstacleCells, nAccelCells, nFluidCells,
		   minRho, maxRho, mass,
		   sqrt ( minU2 ), sqrt ( maxU2 )
	);
}

/*############################################################################*/

static void storeValue ( FILE *file, OUTPUT_PRECISION *v ) {
	const int litteBigEndianTest = 1;
	if  ( (*((unsigned char *)&litteBigEndianTest)) == 0 ) { /* big endian */
		const char *vPtr = (char*) v;
		char buffer[sizeof ( OUTPUT_PRECISION )];

		for ( int i = 0; i < sizeof ( OUTPUT_PRECISION ); ++i ) {
			buffer[i] = vPtr[sizeof ( OUTPUT_PRECISION ) - i - 1];
		}

		fwrite ( buffer, sizeof ( OUTPUT_PRECISION ), 1, file );
	} else { /* little endian */
		fwrite ( v, sizeof ( OUTPUT_PRECISION ), 1, file );
	}
}

/*############################################################################*/

static void loadValue ( FILE *file, OUTPUT_PRECISION *v ) {
	const int litteBigEndianTest = 1;
	if ( (*((unsigned char *)&litteBigEndianTest)) == 0 ) { /* big endian */
		char *vPtr = (char*) v;
		char buffer[sizeof ( OUTPUT_PRECISION )];

		fread ( buffer, sizeof ( OUTPUT_PRECISION ), 1, file );

		for ( int i = 0; i < sizeof ( OUTPUT_PRECISION ); ++i ) {
			vPtr[i] = buffer[sizeof ( OUTPUT_PRECISION ) - i - 1];
		}
	}
	else { /* little endian */
		fread ( v, sizeof ( OUTPUT_PRECISION ), 1, file );
	}
}

/*############################################################################*/

// void LBM_storeVelocityField ( LBM_Grid grid, const char *filename, const int binary ) {
void LBM_storeVelocityField ( float *grid, const char *filename, const int binary ) {
	FILE *file = fopen(filename, (binary ? "wb" : "w"));

	for ( int z = 0; z < SIZE_Z; ++z ) {
		for ( int y = 0; y < SIZE_Y; ++y ) {
			for ( int x = 0; x < SIZE_X; ++x ) {
				float rho = GRID_ENTRY ( grid, x, y, z, C ) + GRID_ENTRY ( grid, x, y, z, N ) + GRID_ENTRY ( grid, x, y, z, S  ) + GRID_ENTRY ( grid, x, y, z, E  ) + GRID_ENTRY ( grid, x, y, z, W  ) + GRID_ENTRY ( grid, x, y, z, T  ) + GRID_ENTRY ( grid, x, y, z, B  ) + GRID_ENTRY ( grid, x, y, z, NE ) + GRID_ENTRY ( grid, x, y, z, NW ) + GRID_ENTRY ( grid, x, y, z, SE ) + GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, NT ) + GRID_ENTRY ( grid, x, y, z, NB ) + GRID_ENTRY ( grid, x, y, z, ST ) + GRID_ENTRY ( grid, x, y, z, SB ) + GRID_ENTRY ( grid, x, y, z, ET ) + GRID_ENTRY ( grid, x, y, z, EB ) + GRID_ENTRY ( grid, x, y, z, WT ) + GRID_ENTRY ( grid, x, y, z, WB );
				float ux  = GRID_ENTRY ( grid, x, y, z, E ) - GRID_ENTRY ( grid, x, y, z, W ) + GRID_ENTRY ( grid, x, y, z, NE ) - GRID_ENTRY ( grid, x, y, z, NW ) + GRID_ENTRY ( grid, x, y, z, SE ) - GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, ET ) + GRID_ENTRY ( grid, x, y, z, EB ) - GRID_ENTRY ( grid, x, y, z, WT ) - GRID_ENTRY ( grid, x, y, z, WB ) ;
				float uy  = GRID_ENTRY ( grid, x, y, z, N ) - GRID_ENTRY ( grid, x, y, z, S ) + GRID_ENTRY ( grid, x, y, z, NE ) + GRID_ENTRY ( grid, x, y, z, NW ) - GRID_ENTRY ( grid, x, y, z, SE ) - GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, NT ) + GRID_ENTRY ( grid, x, y, z, NB ) - GRID_ENTRY ( grid, x, y, z, ST ) - GRID_ENTRY ( grid, x, y, z, SB ) ;
				float uz  = GRID_ENTRY ( grid, x, y, z, T ) - GRID_ENTRY ( grid, x, y, z, B ) + GRID_ENTRY ( grid, x, y, z, NT ) - GRID_ENTRY ( grid, x, y, z, NB ) + GRID_ENTRY ( grid, x, y, z, ST ) - GRID_ENTRY ( grid, x, y, z, SB ) + GRID_ENTRY ( grid, x, y, z, ET ) - GRID_ENTRY ( grid, x, y, z, EB ) + GRID_ENTRY ( grid, x, y, z, WT ) - GRID_ENTRY ( grid, x, y, z, WB ) ;

				ux /= rho;
				uy /= rho;
				uz /= rho;

				if  ( binary ) {
					/*
					fwrite (  &ux, sizeof (  ux  ) , 1, file  ) ;
					fwrite (  &uy, sizeof (  uy  ) , 1, file  ) ;
					fwrite (  &uz, sizeof (  uz  ) , 1, file  ) ;
					*/
					storeValue ( file, &ux ) ;
					storeValue ( file, &uy ) ;
					storeValue ( file, &uz ) ;
				} else {
					fprintf ( file, "%e %e %e\n", ux, uy, uz ) ;
				}
			}
		}
	}

	fclose ( file );
}

/*############################################################################*/

// void LBM_compareVelocityField ( LBM_Grid grid, const char *filename, const int binary ) {
void LBM_compareVelocityField ( float *grid, const char *filename, const int binary ) {
	float maxDiff2 = -1e+30;

	FILE *file = fopen ( filename, ( binary ? "rb" : "r" ) );

	for ( int z = 0; z < SIZE_Z; ++z ) {
		for ( int y = 0; y < SIZE_Y; ++y ) {
			for ( int x = 0; x < SIZE_X; ++x ) {
				float rho = GRID_ENTRY ( grid, x, y, z, C ) + GRID_ENTRY ( grid, x, y, z, N ) + GRID_ENTRY ( grid, x, y, z, S  ) + GRID_ENTRY ( grid, x, y, z, E  ) + GRID_ENTRY ( grid, x, y, z, W  ) + GRID_ENTRY ( grid, x, y, z, T  ) + GRID_ENTRY ( grid, x, y, z, B  ) + GRID_ENTRY ( grid, x, y, z, NE ) + GRID_ENTRY ( grid, x, y, z, NW ) + GRID_ENTRY ( grid, x, y, z, SE ) + GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, NT ) + GRID_ENTRY ( grid, x, y, z, NB ) + GRID_ENTRY ( grid, x, y, z, ST ) + GRID_ENTRY ( grid, x, y, z, SB ) + GRID_ENTRY ( grid, x, y, z, ET ) + GRID_ENTRY ( grid, x, y, z, EB ) + GRID_ENTRY ( grid, x, y, z, WT ) + GRID_ENTRY ( grid, x, y, z, WB );
				float ux  = GRID_ENTRY ( grid, x, y, z, E ) - GRID_ENTRY ( grid, x, y, z, W ) + GRID_ENTRY ( grid, x, y, z, NE ) - GRID_ENTRY ( grid, x, y, z, NW ) + GRID_ENTRY ( grid, x, y, z, SE ) - GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, ET ) + GRID_ENTRY ( grid, x, y, z, EB ) - GRID_ENTRY ( grid, x, y, z, WT ) - GRID_ENTRY ( grid, x, y, z, WB );
				float uy  = GRID_ENTRY ( grid, x, y, z, N ) - GRID_ENTRY ( grid, x, y, z, S ) + GRID_ENTRY ( grid, x, y, z, NE ) + GRID_ENTRY ( grid, x, y, z, NW ) - GRID_ENTRY ( grid, x, y, z, SE ) - GRID_ENTRY ( grid, x, y, z, SW ) + GRID_ENTRY ( grid, x, y, z, NT ) + GRID_ENTRY ( grid, x, y, z, NB ) - GRID_ENTRY ( grid, x, y, z, ST ) - GRID_ENTRY ( grid, x, y, z, SB );
				float uz  = GRID_ENTRY ( grid, x, y, z, T ) - GRID_ENTRY ( grid, x, y, z, B ) + GRID_ENTRY ( grid, x, y, z, NT ) - GRID_ENTRY ( grid, x, y, z, NB ) + GRID_ENTRY ( grid, x, y, z, ST ) - GRID_ENTRY ( grid, x, y, z, SB ) + GRID_ENTRY ( grid, x, y, z, ET ) - GRID_ENTRY ( grid, x, y, z, EB ) + GRID_ENTRY ( grid, x, y, z, WT ) - GRID_ENTRY ( grid, x, y, z, WB );
				ux /= rho;
				uy /= rho;
				uz /= rho;

				OUTPUT_PRECISION fileUx, fileUy, fileUz;
				if ( binary ) {
					loadValue ( file, &fileUx );
					loadValue ( file, &fileUy );
					loadValue ( file, &fileUz );
				} else {
					if ( sizeof ( OUTPUT_PRECISION ) == sizeof ( double ) ) {
						fscanf ( file, "%lf %lf %lf\n", &fileUx, &fileUy, &fileUz );
					} else {
						fscanf ( file, "%f %f %f\n", &fileUx, &fileUy, &fileUz );
					}
				}

				float dUx = ux - fileUx;
				float dUy = uy - fileUy;
				float dUz = uz - fileUz;
				float diff2 = dUx * dUx + dUy * dUy + dUz * dUz;
				if (diff2 > maxDiff2) {
					maxDiff2 = diff2;
				}
			}
		}
	}

#if defined ( SPEC_CPU )
	printf ( "LBM_compareVelocityField: maxDiff = %e  \n\n", sqrt ( maxDiff2 ) );
#else
	printf ( "LBM_compareVelocityField: maxDiff = %e  ==>  %s\n\n", sqrt ( maxDiff2 ), sqrt ( maxDiff2 ) > 1e-5 ? "##### ERROR #####" : "OK" );
#endif
	fclose ( file );
}
