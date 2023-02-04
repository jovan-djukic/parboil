/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "model.h"

#include <stdlib.h>

#include "args.h"

int doCompute ( struct cartesian *data1, int n1, struct cartesian *data2, int n2, int doSelf, long long *data_bins, int nbins, float *binb ) {
	if ( doSelf ) {
		n2    = n1;
		data2 = data1;
	}

	const unsigned int iLimit = doSelf ? ( n1 - 1 ) : n1;
	#pragma acc parallel loop present(data1[0:n1],data2[0:n2],binb[0:(nbins + 1)], data_bins[0:(nbins + 2)]) 
	for ( unsigned int i = 0; i < iLimit; ++i ) {
		const unsigned int jInitialValue = doSelf ? ( i + 1 ) : 0;
		for ( unsigned int j = jInitialValue; j < n2; ++j ) {
			const float xi = data1[i].x;
			const float yi = data1[i].y;
			const float zi = data1[i].z;
			const float dot = xi * data2[j].x + yi * data2[j].y + zi * data2[j].z;

			// run binary search
			int min = 0;
			int max = nbins;
			int index;

			while ( max > ( min + 1 ) ) {
				int k = ( min + max ) / 2;
				if ( dot >= binb[k] ) {
					max = k;
				} else {
					min = k;
				}
			};

			if ( dot >= binb[min] ) {
				// #pragma acc atomic update	
				// data_bins[min] += 1; /*k = min;*/
				index = min;
			} else if ( dot < binb[max] ) {
				// #pragma acc atomic update	
				// data_bins[max + 1] += 1; /*k = max+1;*/
				index = max + 1;
			} else {
				// #pragma acc atomic update	
				// data_bins[max] += 1; /*k = max;*/
				index = max;
			}

			#pragma acc atomic update
			data_bins[index] += 1;
		}
	}

	return 0;
}

void doComputeWrapper ( options args, struct cartesian *data, struct cartesian *random, long long *DD, long long *RRS, long long *DRS, int nbins, float *binb ) {
    #pragma acc parallel loop collapse(3) copyin(data[0:args.npoints], random[0:(args.npoints * args.random_count)], binb[0:(nbins + 1)], args) copy(DD[0:(nbins + 2)], RRS[0:( ( nbins + 2 ) * args.random_count )], DRS[0:( ( nbins + 2 ) * args.random_count )])
    for ( int rf = 0; rf < ( 2 * args.random_count + 1 ); ++rf ) {
		for ( unsigned int i = 0; i < args.npoints; ++i ) {
			for ( unsigned int j = 0; j < args.npoints; ++j ) {
				// rf == 0 => DD
				// rf is odd => RRS
				// else its DRS
				int doSelf = ( rf == 0 ) || ( ( rf & 1 ) != 0 );

				int skip = ( doSelf == 1 ) && ( ( i == ( args.npoints - 1 ) ) || ( j < ( i + 1 ) ) ); 

				if ( skip == 0 ) {
					struct cartesian *data1 = rf == 0 ? data : ( ( rf & 1 ) != 0 ? ( random + ( rf / 2 ) * args.npoints ) : data );				
					struct cartesian *data2 = rf == 0 ? data : ( ( rf & 1 ) != 0 ? ( random + ( rf / 2 ) * args.npoints ) : ( random + ( ( rf - 1 ) / 2 ) * args.npoints ) );				

					long long *data_bins = rf == 0 ? DD : ( ( rf & 1 ) != 0 ? ( RRS + ( rf / 2 ) * ( nbins + 2 ) ) : ( DRS + ( ( rf - 1 ) / 2 ) * ( nbins + 2 ) ) );

					const float xi = data1[i].x;
					const float yi = data1[i].y;
					const float zi = data1[i].z;
					const float dot = xi * data2[j].x + yi * data2[j].y + zi * data2[j].z;

					// run binary search
					int min = 0;
					int max = nbins;
					int index;

					while ( max > ( min + 1 ) ) {
						int k = ( min + max ) / 2;
						if ( dot >= binb[k] ) {
							max = k;
						} else {
							min = k;
						}
					};

					if ( dot >= binb[min] ) {
						index = min;
					} else if ( dot < binb[max] ) {
						index = max + 1;
					} else {
						index = max;
					}

					#pragma acc atomic update
					data_bins[index] += 1;
				}
			}
		}
	}

	// reduce RRS and DRS
	for ( int i = 0; i < ( nbins + 2 ); ++i ) {
		for ( int j = 1; j < args.random_count; ++j ) {
			RRS[i] += RRS[i + j * ( nbins + 2 )];
			DRS[i] += DRS[i + j * ( nbins + 2 )];
		}
	}
}