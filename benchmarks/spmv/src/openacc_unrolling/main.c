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

#include "file.h"
#include "convert_dataset.h"

int main ( int argc, char **argv ) {
	struct pb_TimerSet timers;
	struct pb_Parameters *parameters;

	printf ( "OPENACC-based sparse matrix vector multiplication****\n" );
	parameters = pb_ReadParameters ( &argc, argv );
	if ( ( parameters->inpFiles[0] == NULL ) || ( parameters->inpFiles[1] == NULL ) ) {
		fprintf ( stderr, "Expecting two input filenames\n" );
		exit ( -1 );
	}

	pb_InitializeTimerSet ( &timers );
	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	//parameters declaration
	int len;
	int depth;
	int dim;
	int pad = 1;
	int nzcnt_len;

	//host memory allocation
	//matrix
	float *h_data;
	int *h_indices;
	int *h_ptr;
	int *h_perm;
	int *h_nzcnt;
	//vector
	float *h_Ax_vector;
	float *h_x_vector;

	//load matrix from files
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );
	//inputData(parameters->inpFiles[0], &len, &depth, &dim,&nzcnt_len,&pad,
	//    &h_data, &h_indices, &h_ptr,
	//    &h_perm, &h_nzcnt);
	int col_count;
	coo_to_jds (
		parameters->inpFiles[0], // bcsstk32.mtx, fidapm05.mtx, jgl009.mtx
		1,						 // row padding
		pad,					 // warp size
		1,						 // pack size
		1,						 // is mirrored?
		0,						 // binary matrix
		1,						 // debug level [0:2]
		&h_data, &h_ptr, &h_nzcnt, &h_indices, &h_perm,
		&col_count, &dim, &len, &nzcnt_len, &depth
	);

	h_Ax_vector = ( float* )malloc ( sizeof ( float ) * dim );
	h_x_vector = ( float* )malloc ( sizeof ( float ) * dim );
	input_vec ( parameters->inpFiles[1], h_x_vector, dim );

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	// main execution
	#pragma acc data copyin(h_nzcnt[0:nzcnt_len], h_ptr[0:depth], h_indices[0:len], h_data[0:len], h_x_vector[0:dim], h_perm[0:dim]) copyout(h_Ax_vector[0:dim])
	{
		for ( unsigned int timeStep = 0; timeStep < 50; ++timeStep ) {
			#pragma acc parallel loop
			for ( unsigned int ix = 0; ix < dim; ++ix ) {
				float sum = 0.0f;
				int	bound = h_nzcnt[ix];

				//prefetch 0
				int ptr = h_ptr[0] + ix;  
				float d = h_data[ptr]; 
				int index = h_indices[ptr];  
				float t = h_x_vector[index];
		
				if ( bound > 1 ) {  //bound >= 2
					//prefetch 1
					ptr = h_ptr[1] + ix;    
					index = h_indices[ptr];  

					float dn = 0.0f;
					float tn = 0.0f;

					for ( int k = 2; k < bound; ++k ) {	
						//prefetch k-1
						dn = h_data[ptr]; 

						//prefetch k
						ptr = h_ptr[k] + ix;    

						int in = h_indices[ptr]; 

						//prefetch k-1
						tn = h_x_vector[index];
						
						//compute k-2
						sum += d * t; 

						//sweep to k
						index = in;  

						//sweep to k-1
						d = dn;
						t = tn; 
					}	
		
					//fetch last
					dn = h_data[ptr];
					tn = h_x_vector[index];
	
					//compute last-1
					sum += d * t; 

					//sweep to last
					d = dn;
					t = tn;
				}
				//compute last
				sum += d * t;  

				h_Ax_vector[h_perm[ix]] = sum; 
			}
		}
	}

	if ( parameters->outFile ) {
		pb_SwitchToTimer ( &timers, pb_TimerID_IO );
		outputData ( parameters->outFile, h_Ax_vector, dim );
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	free ( h_data );
	free ( h_indices );
	free ( h_ptr );
	free ( h_perm );
	free ( h_nzcnt );
	free ( h_Ax_vector );
	free ( h_x_vector );
	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );

	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );

	return 0;
}
