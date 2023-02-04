/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <parboil.h>

#include "sad.h"
#include "file.h"
#include "image.h"

static unsigned short *load_sads ( char *filename );
static void write_sads ( char *filename, int image_width_macroblocks, int image_height_macroblocks, unsigned short *sads );
static void write_sads_directly ( char *filename, int width, int height, unsigned short *sads );

/* FILE I/O */
unsigned short* load_sads ( char *filename ) {

	FILE *file = fopen ( filename, "r" );

	if ( !file ) {
		fprintf ( stderr, "Cannot find file '%s'\n", filename );
		exit ( -1 );
	}

	/* Read image dimensions (measured in macroblocks) */
	int width = read16u ( file );
	int height = read16u ( file );

	/* Read SAD values.  Only interested in the 4x4 SAD values, which are
	 * at the end of the file. */
	int sads_per_block = MAX_POS_PADDED * ( width * height );
	fseek ( file, 25 * sads_per_block * sizeof ( unsigned short ), SEEK_CUR );

	unsigned short *sads = ( unsigned short* ) malloc ( sads_per_block * 16 * sizeof ( unsigned short ) );
	fread ( sads, sizeof ( unsigned short ), sads_per_block * 16, file );
	fclose ( file );

	return sads;
}

/* Compare the reference SADs to the expected SADs. */
void check_sads ( unsigned short *sads_reference, unsigned short *sads_computed, int image_size_macroblocks ) {
	/* Check the 4x4 SAD values.  These are in sads_reference.
	 * Ignore the data at the beginning of sads_computed. */
	sads_computed += 25 * MAX_POS_PADDED * image_size_macroblocks;

	for ( int block = 0; block < image_size_macroblocks; ++block ) {
		for ( int subblock = 0; subblock < 16; ++subblock ) {
			for ( int sad_index = 0; sad_index < MAX_POS; ++sad_index ) {
				int index = ( block * 16 + subblock ) * MAX_POS_PADDED + sad_index;

				if ( sads_reference[index] != sads_computed[index] ) {
#if 0
		  /* Print exactly where the mismatch was seen */
		  printf("M %3d %2d %4d (%d = %d)\n", block, subblock, sad_index, sads_reference[index], sads_computed[index]);
#else
					goto mismatch;
#endif
				}
			}
		}
	}

	printf ( "Success.\n" );
	return;

mismatch:
	printf ( "Computed SADs do not match expected values.\n" );
}

/* Extract the SAD data for a particular block type for a particular
 * macroblock from the array of SADs of that block type. */
static inline void write_subblocks ( FILE *outfile, unsigned short *subblock_array, int macroblock, int count ) {
	for ( int block = 0; block < count; ++block ) {
		unsigned short *vector = subblock_array + ( block + macroblock * count ) * MAX_POS_PADDED;

		/* Write all SADs for this sub-block */
		for ( int position = 0; position < MAX_POS; ++position ) {
			write16u ( outfile, *vector++ );
		}
	}
}

/* Write some SAD data to a file for output checking.
 *
 * All SAD values for six rows of macroblocks are written.
 * The six rows consist of the top two, middle two, and bottom two image rows.
 */
void write_sads ( char *filename, int mb_width, int mb_height, unsigned short *sads ) {
	FILE *file = fopen ( filename, "w" );
	int mbs = mb_width * mb_height;
	int row_indices[6] = {
		0, 1,
		 mb_height / 2 - 1, mb_height / 2,
		 mb_height - 2, mb_height - 1
	};

	if ( !file ) {
		fprintf ( stderr, "Cannot open output file\n" );
		exit ( -1 );
	}

	/* Write the number of output macroblocks */
	write32u ( file, mb_width * 6 );

	/* Write zeros */
	write32u ( file, 0 );

	/* Each row */
	for ( int row_index = 0; row_index < 6; ++row_index ) {
		int row = row_indices[row_index];

		/* Each block in row */
		for ( int block = mb_width * row; block < mb_width * ( row + 1 ); ++block ) {
			/* Write SADs for all sub-block types */
			for ( int blocktype = 1; blocktype <= 7; ++blocktype ) {
				write_subblocks ( 
					file,
					sads + SAD_TYPE_IX ( blocktype, mbs ),
					block,
					SAD_TYPE_CT ( blocktype )
				);
			}
		}
	}

	fclose ( file );
}

/* FILE I/O for debugging */
static void write_sads_directly ( char *filename, int width, int height, unsigned short *sads ) {
	FILE *file = fopen ( filename, "w" );

	write16u ( file, width );
	write16u ( file, height );

	for ( int n = 0; n < 41 * MAX_POS_PADDED * ( width * height ); ++n ) {
		write16u ( file, sads[n] );
	}

	fclose ( file );
}

static void print_test_sad_vector ( unsigned short *base, int macroblock, int count ) {
	const int searchpos = 17 * 33 + 17;
	for ( int n = 0; n < count; ++n ) {
		printf ( " %d", base[( count * macroblock + n ) * MAX_POS_PADDED + searchpos] );
	}
}

static void print_test_sads ( unsigned short *sads_computed, int mbs ) {
	int macroblock = 5;

	for ( int blocktype = 1; blocktype <= 7; ++blocktype ) {
		printf ( "%d:", blocktype );
		print_test_sad_vector ( sads_computed + SAD_TYPE_IX ( blocktype, mbs ), macroblock, SAD_TYPE_CT ( blocktype ) );
		puts ( "\n" );
	}
}

/* MAIN */
int main ( int argc, char **argv ) {
	struct pb_TimerSet timers;
	pb_InitializeTimerSet(&timers);
	struct pb_Parameters *parameters = pb_ReadParameters(&argc, argv);

	if ( pb_Parameters_CountInputs ( parameters ) != 2 ) {
		fprintf ( stderr, "Expecting two input filenames\n" );
		exit ( -1 );
	}

	/* Read input files */
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	struct image_i16 *ref_image = load_image ( parameters->inpFiles[0] );
	struct image_i16 *cur_image = load_image ( parameters->inpFiles[1] );

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	if ( ( ref_image->width != cur_image->width ) || ( ref_image->height != cur_image->height ) ) {
		fprintf ( stderr, "Input images must be the same size\n" );
		exit ( -1 );
	}

	if ( ( ref_image->width % 16 ) || ( ref_image->height % 16 ) ) {
		fprintf ( stderr, "Input image size must be an integral multiple of 16\n" );
		exit ( -1 );
	}

	/* Compute parameters, allocate memory */
	int image_size_bytes         = ref_image->width * ref_image->height * sizeof ( short );
	int image_width_macroblocks  = ref_image->width / 16;
	int image_height_macroblocks = ref_image->height / 16;
	int image_size_macroblocks   = image_width_macroblocks * image_height_macroblocks;

	/* SADs generated by the program */
	unsigned short *sads_computed = ( unsigned short* ) malloc ( 41 * MAX_POS_PADDED * image_size_macroblocks * sizeof ( short ) );

	/* Run the kernel code */
	sad4_cpu ( 
		sads_computed,
		( unsigned short* ) cur_image->data,
		( unsigned short* ) ref_image->data,
		ref_image->width / 16, ref_image->height / 16
	);

	larger_sads ( sads_computed, image_size_macroblocks );

	/* Print output */
	if ( parameters->outFile ) {
		pb_SwitchToTimer ( &timers, pb_TimerID_IO );

		write_sads (
			parameters->outFile,
			image_width_macroblocks,
			image_height_macroblocks,
			sads_computed
		);

		pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
	}

#if 0 /* Debugging */
  print_test_sads(sads_computed, image_size_macroblocks);
  write_sads_directly("sad-debug.bin",
		      ref_image->width / 16, ref_image->height / 16,
		      sads_computed);
#endif

	/* Free memory */
	free ( sads_computed );
	free_image ( ref_image );
	free_image ( cur_image );

	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );
	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );

	return 0;
}
