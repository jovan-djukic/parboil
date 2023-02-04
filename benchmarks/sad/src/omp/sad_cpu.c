/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdlib.h>
#include "sad.h"

void sad4_cpu ( unsigned short *blk_sad, unsigned short *frame, unsigned short *ref, int mb_width, int mb_height, int image_size_macroblocks ) {
	/* For each block */
	#pragma omp parallel for
	for ( int mb_y = 0; mb_y < mb_height; ++mb_y ) {
		for ( int mb_x = 0; mb_x < mb_width; ++mb_x ) {
		  	/* Each search position */
		  	for ( int position_y = -SEARCH_RANGE; position_y <= SEARCH_RANGE; ++position_y ) {
		    	for ( int position_x = -SEARCH_RANGE; position_x <= SEARCH_RANGE; ++position_x ) {
		      		/* Each 4x4 block in the macroblock */
		      		for ( int block_y = 0; block_y < 4; ++block_y ) {
						for ( int block_x = 0; block_x < 4; ++block_x ) {
							/* Go to the starting offset in blk_sad */
							unsigned short *macroblock_sad =  blk_sad + SAD_TYPE_7_IX ( mb_width * mb_height ) + ( mb_y * mb_width + mb_x ) * ( SAD_TYPE_7_CT * MAX_POS_PADDED );
							unsigned short *current_frame = frame + mb_y * 256 * mb_width + mb_x * 16;

							int frame_y	= mb_y * 16;
							int	frame_x = mb_x * 16;

							int width  = mb_width * 16;
							int height = mb_height * 16;

		 					/* search position */
		 					int position = ( position_y + SEARCH_RANGE ) * ( 2 * SEARCH_RANGE + 1 ) + ( position_x + SEARCH_RANGE );

			  				unsigned short sad = 0;

			  				/* Each pixel */
			  				for ( int y = 0; y < 4; ++y ) {
			    				for ( int x = 0; x < 4; ++x ) {
			    				  /* Get reference pixel coordinate, clipped to image boundary */
									int reference_x = frame_x + position_x + ( block_x * 4 ) + x;
									int reference_y = frame_y + position_y + ( block_y * 4 ) + y;

			    				  	if ( reference_x < 0 ) {
								  		reference_x = 0;
								  	} else if ( reference_x >= width ) {
								  		reference_x = width - 1;
								  	}

			    					if ( reference_y < 0 ) {
										reference_y = 0;
									} else if ( reference_y >= height ) {
										reference_y = height - 1;
									}

									unsigned int b = ref[reference_y * width + reference_x];
									unsigned int a = current_frame[( block_y * 4 + y ) * width + ( block_x * 4 + x )];

			    				  	sad += abs ( a - b );
			    				}
							}

							/* Save the SAD */
							macroblock_sad[MAX_POS_PADDED*( 4 * block_y + block_x ) + position] = sad;
						}
		      		}
		    	}
			}
		}
	}
}

void add_vectors ( unsigned short *x, unsigned short *y, unsigned short *z ) {
	for ( int count = 0; count < (MAX_POS + 1) / 2; ++count, z += 2, x += 2, y += 2) {
		unsigned int *x_cast = ( unsigned int* ) x;
		unsigned int *y_cast = ( unsigned int* ) y;
		unsigned int *z_cast = ( unsigned int* ) z;		

		*z_cast = *x_cast + *y_cast;
	}
}

void larger_sads ( unsigned short *sads, int mbs ) {
	/* x and y are inputs to vector addition */
	/* z is the output of vector addition */

	const int size = 41 * MAX_POS_PADDED * mbs;
	#pragma acc parallel loop present(sads[0:size])
	for ( int macroblock = 0; macroblock < mbs; ++macroblock ) {
		/* Block type 6 */
		for ( int block_y = 0; block_y < 2; ++block_y ) {
			for ( int block_x = 0; block_x < 4; ++block_x ) {
				unsigned short *x = sads + SAD_TYPE_7_IX ( mbs ) + macroblock * SAD_TYPE_7_CT * MAX_POS_PADDED + ( 8 * block_y + block_x ) * MAX_POS_PADDED;
				unsigned short *y = x + 4 * MAX_POS_PADDED;
				unsigned short *z = sads + SAD_TYPE_6_IX ( mbs ) + macroblock * SAD_TYPE_6_CT * MAX_POS_PADDED + ( 4 * block_y + block_x ) * MAX_POS_PADDED;

				add_vectors ( x, y, z );
			}
		}

		/* Block type 5 */
		for ( int block_y = 0; block_y < 4; ++block_y ) {
			for ( int block_x = 0; block_x < 2; ++block_x ) {
				unsigned short *x = sads + SAD_TYPE_7_IX ( mbs ) + macroblock * SAD_TYPE_7_CT * MAX_POS_PADDED + ( 4 * block_y + 2 * block_x ) * MAX_POS_PADDED;
				unsigned short *y = x + MAX_POS_PADDED;
				unsigned short *z = sads + SAD_TYPE_5_IX ( mbs ) + macroblock * SAD_TYPE_6_CT * MAX_POS_PADDED + ( 2 * block_y + block_x ) * MAX_POS_PADDED;

				add_vectors ( x, y, z );
			}
		}

		/* Block type 4 */
		for ( int block_y = 0; block_y < 2; ++block_y ) {
			for ( int block_x = 0; block_x < 2; ++block_x ) {
				unsigned short *x = sads + SAD_TYPE_5_IX ( mbs ) + macroblock * SAD_TYPE_5_CT * MAX_POS_PADDED + ( 4 * block_y + block_x ) * MAX_POS_PADDED;
				unsigned short *y = x + 2 * MAX_POS_PADDED;
				unsigned short *z = sads + SAD_TYPE_4_IX ( mbs ) + macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED + ( 2 * block_y + block_x ) * MAX_POS_PADDED;

				add_vectors ( x, y, z );
			}
		}

		/* Block type 3 */
		{
			unsigned short *x = sads + SAD_TYPE_4_IX ( mbs ) + macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED;
			unsigned short *y = x + 2 * MAX_POS_PADDED;
			unsigned short *z = sads + SAD_TYPE_3_IX ( mbs ) + macroblock * SAD_TYPE_3_CT * MAX_POS_PADDED;

			add_vectors ( x, y, z );
		}

		{
			unsigned short *x = sads + SAD_TYPE_4_IX ( mbs ) + macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED + MAX_POS_PADDED;
			unsigned short *y = x + 2 * MAX_POS_PADDED;
			unsigned short *z = sads + SAD_TYPE_3_IX ( mbs ) + macroblock * SAD_TYPE_3_CT * MAX_POS_PADDED + MAX_POS_PADDED;

			add_vectors ( x, y, z );
		}

		/* Block type 2 */
		{
			unsigned short *x = sads + SAD_TYPE_4_IX ( mbs ) + macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED;
			unsigned short *y = x + MAX_POS_PADDED;
			unsigned short *z = sads + SAD_TYPE_2_IX ( mbs ) + macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED;

			add_vectors ( x, y, z );
		}

		{
			unsigned short *x = sads + SAD_TYPE_4_IX ( mbs ) + macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED + 2 * MAX_POS_PADDED;
			unsigned short *y = x + MAX_POS_PADDED;
			unsigned short *z = sads + SAD_TYPE_2_IX ( mbs ) + macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED + MAX_POS_PADDED;

			add_vectors ( x, y, z );
		}

		/* Block type 1 */
		{
			unsigned short *x = sads + SAD_TYPE_2_IX ( mbs ) + macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED;
			unsigned short *y = x + MAX_POS_PADDED;
			unsigned short *z = sads + SAD_TYPE_1_IX ( mbs ) + macroblock * SAD_TYPE_1_CT * MAX_POS_PADDED;

			add_vectors ( x, y, z );
		}
	}
}
