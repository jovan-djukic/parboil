/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/* The optimized CPU-only version of the code, for performance measurement. */

#include <stdlib.h>
#include "sad.h"

static short line[16];

static void sad4_one_macroblock(
	unsigned short *blk_sad,
	unsigned short *frame,
	unsigned short *ref,
	int org_y,
	int org_x,
	int mb_width,
	int mb_height
	);

static short *create_padded_row ( short *ref, int y, int x, int height, int width );

void sad4_cpu(unsigned short *blk_sad, unsigned short *frame, unsigned short *ref, int mb_width, int mb_height) {
	/* Go to the starting offset in blk_sad */
	blk_sad += SAD_TYPE_7_IX(mb_width * mb_height);

	/* For each block */
	for ( int mb_y = 0, frame_yoff = 0; mb_y < mb_height; ++mb_y, frame_yoff += 256 * mb_width ) {
		for ( int mb_x = 0; mb_x < mb_width; ++mb_x ) {
			sad4_one_macroblock(
				blk_sad + ( mb_y * mb_width + mb_x ) * ( SAD_TYPE_7_CT * MAX_POS_PADDED ),
				frame + frame_yoff + mb_x * 16,
				ref,
				mb_y * 16,
				mb_x * 16,
				mb_width,
				mb_height
			);
		}
	}
}

void sad4_one_macroblock(
	unsigned short *macroblock_sad,
	unsigned short *frame,
	unsigned short *ref,
	int frame_y,
	int frame_x,
	int mb_width,
	int mb_height
) {
	/* current macroblock in frame */
	unsigned short frame_mb[256] = { 0 };

	int width  = mb_width * 16;
	int height = mb_height * 16;

	int max_width  = width - 17;
	int max_height = height - 17;


	/* Make a local copy of frame */
	for ( int y = 0; y < 16; ++y ) {
		for ( int x = 0; x < 16; ++x ) {
			frame_mb[16 * y + x] = frame[width * y + x];
		}
	}

	int range_partly_outside = 0;
	int frame_x_in_search_range = ( frame_x >= SEARCH_RANGE ) && ( frame_x <= width - SEARCH_RANGE - 17 );
	int frame_y_in_search_range = ( frame_y >= SEARCH_RANGE ) && ( frame_y <= height - SEARCH_RANGE - 17 );
	if ( frame_x_in_search_range && frame_y_in_search_range ) {
		range_partly_outside = 0;
	} else {
		range_partly_outside = 1;
	}

	/* search position */
	int position = 0;

	/* Each search position */
	for ( int position_y = -SEARCH_RANGE; position_y <= SEARCH_RANGE; ++position_y ) {
		for ( int position_x = -SEARCH_RANGE; position_x <= SEARCH_RANGE; ++position_x, ++position ) {
			short *current_pointer = frame_mb;

			unsigned short *sad_line = macroblock_sad;

			int absolute_y = frame_y + position_y;
			int absolute_x = frame_x + position_x;

			int absolute_y_in_range = absolute_y >= 0 && absolute_y <= max_height;
			int absolute_x_in_range = absolute_x >= 0 && absolute_x <= max_width;
			int do_bounds_check = range_partly_outside && !( absolute_y_in_range && absolute_x_in_range );

			for ( int block_y = 0; block_y < 4; ++block_y ) {
				int sad0 = 0;
				int sad1 = 0;
				int sad2 = 0;
				int sad3 = 0;

				for ( int y = 0; y < 4; ++y ) {
					short *reference_pointer = 	do_bounds_check
								 				? create_padded_row ( ref, absolute_y, absolute_x, height, width )
								 				: ( short* ) ( ref + absolute_y * width + absolute_x );

					absolute_y++;

					sad0 += abs ( *reference_pointer++ - *current_pointer++ );
					sad0 += abs ( *reference_pointer++ - *current_pointer++ );
					sad0 += abs ( *reference_pointer++ - *current_pointer++ );
					sad0 += abs ( *reference_pointer++ - *current_pointer++ );
					sad1 += abs ( *reference_pointer++ - *current_pointer++ );
					sad1 += abs ( *reference_pointer++ - *current_pointer++ );
					sad1 += abs ( *reference_pointer++ - *current_pointer++ );
					sad1 += abs ( *reference_pointer++ - *current_pointer++ );
					sad2 += abs ( *reference_pointer++ - *current_pointer++ );
					sad2 += abs ( *reference_pointer++ - *current_pointer++ );
					sad2 += abs ( *reference_pointer++ - *current_pointer++ );
					sad2 += abs ( *reference_pointer++ - *current_pointer++ );
					sad3 += abs ( *reference_pointer++ - *current_pointer++ );
					sad3 += abs ( *reference_pointer++ - *current_pointer++ );
					sad3 += abs ( *reference_pointer++ - *current_pointer++ );
					sad3 += abs ( *reference_pointer++ - *current_pointer++ );
				}

				sad_line[                     position] = sad0;
				sad_line[MAX_POS_PADDED     + position] = sad1;
				sad_line[MAX_POS_PADDED * 2 + position] = sad2;
				sad_line[MAX_POS_PADDED * 3 + position] = sad3;

				sad_line += MAX_POS_PADDED * 4;
			}
		}
	}
}

/* Return a row of 16 pixels starting at offset (x, y).  The row may lie
 * partly outside the image, in which case an appropriate row will be
 * constructed in 'line' and returned.  Otherwise, a reference to the
 * image is returned. */
static short * create_padded_row ( short *ref, int y, int x, int height, int width ) {
	if ( y < 0 ) {
		y = 0;
	} else if ( y >= height ) {
		y = height - 1;
	}

	short *row = ref + y * width;

	if ( ( x >= 0 ) && ( x <= width - 16 ) ) {
		return row + x;
	}

	int i = 0;
	/* Pad left edge of image */
	for ( ; ( x < 0 ) && ( i < 16 ); ++x, ++i ) {
		line[i] = row[0];
	}

	/* Copy row from image */
	for ( ; ( x < width ) && ( i < 16 ); ++x, ++i ) {
		line[i] = row[x];
	}

	/* Pad right edge of image */
	for ( ; i < 16; ++x, ++i ) {
		line[i] = row[width - 1];
	}

	return line;
}

void larger_sads(unsigned short *sads, int mbs)
{
	int macroblock;
	int block_x, block_y;
	unsigned short *x, *y; /* inputs to vector addition */
	unsigned short *z;	   /* output of vector addition */
	int count;

#if 0
#define ADD_VECTORS()                         \
	for (count = 0; count < MAX_POS; count++) \
	*z++ = *x++ + *y++
#else
#define ADD_VECTORS()                                                           \
	for (count = 0; count < (MAX_POS + 1) / 2; count++, z += 2, x += 2, y += 2) \
	*(unsigned int *)z = *(unsigned int *)x + *(unsigned int *)y
#endif

	for (macroblock = 0; macroblock < mbs; macroblock++)
	{
		/* Block type 6 */
		for (block_y = 0; block_y < 2; block_y++)
			for (block_x = 0; block_x < 4; block_x++)
			{
				x = sads + SAD_TYPE_7_IX(mbs) +
					macroblock * SAD_TYPE_7_CT * MAX_POS_PADDED +
					(8 * block_y + block_x) * MAX_POS_PADDED;
				y = x + 4 * MAX_POS_PADDED;
				z = sads + SAD_TYPE_6_IX(mbs) +
					macroblock * SAD_TYPE_6_CT * MAX_POS_PADDED +
					(4 * block_y + block_x) * MAX_POS_PADDED;

				ADD_VECTORS();
			}

		/* Block type 5 */
		for (block_y = 0; block_y < 4; block_y++)
			for (block_x = 0; block_x < 2; block_x++)
			{
				x = sads + SAD_TYPE_7_IX(mbs) +
					macroblock * SAD_TYPE_7_CT * MAX_POS_PADDED +
					(4 * block_y + 2 * block_x) * MAX_POS_PADDED;
				y = x + MAX_POS_PADDED;
				z = sads + SAD_TYPE_5_IX(mbs) +
					macroblock * SAD_TYPE_6_CT * MAX_POS_PADDED +
					(2 * block_y + block_x) * MAX_POS_PADDED;

				ADD_VECTORS();
			}

		/* Block type 4 */
		for (block_y = 0; block_y < 2; block_y++)
			for (block_x = 0; block_x < 2; block_x++)
			{
				x = sads + SAD_TYPE_5_IX(mbs) +
					macroblock * SAD_TYPE_5_CT * MAX_POS_PADDED +
					(4 * block_y + block_x) * MAX_POS_PADDED;
				y = x + 2 * MAX_POS_PADDED;
				z = sads + SAD_TYPE_4_IX(mbs) +
					macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED +
					(2 * block_y + block_x) * MAX_POS_PADDED;

				ADD_VECTORS();
			}

		/* Block type 3 */
		x = sads + SAD_TYPE_4_IX(mbs) +
			macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED;
		y = x + 2 * MAX_POS_PADDED;
		z = sads + SAD_TYPE_3_IX(mbs) +
			macroblock * SAD_TYPE_3_CT * MAX_POS_PADDED;

		ADD_VECTORS();

		x = sads + SAD_TYPE_4_IX(mbs) +
			macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED +
			MAX_POS_PADDED;
		y = x + 2 * MAX_POS_PADDED;
		z = sads + SAD_TYPE_3_IX(mbs) +
			macroblock * SAD_TYPE_3_CT * MAX_POS_PADDED +
			MAX_POS_PADDED;

		ADD_VECTORS();

		/* Block type 2 */
		x = sads + SAD_TYPE_4_IX(mbs) +
			macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED;
		y = x + MAX_POS_PADDED;
		z = sads + SAD_TYPE_2_IX(mbs) +
			macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED;

		ADD_VECTORS();

		x = sads + SAD_TYPE_4_IX(mbs) +
			macroblock * SAD_TYPE_4_CT * MAX_POS_PADDED +
			2 * MAX_POS_PADDED;
		y = x + MAX_POS_PADDED;
		z = sads + SAD_TYPE_2_IX(mbs) +
			macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED +
			MAX_POS_PADDED;

		ADD_VECTORS();

		/* Block type 1 */
		x = sads + SAD_TYPE_2_IX(mbs) +
			macroblock * SAD_TYPE_2_CT * MAX_POS_PADDED;
		y = x + MAX_POS_PADDED;
		z = sads + SAD_TYPE_1_IX(mbs) +
			macroblock * SAD_TYPE_1_CT * MAX_POS_PADDED;

		ADD_VECTORS();
	}
}
