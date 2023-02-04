#include <parboil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

#include <cuda.h>

#define MAX_THREADS_PER_BLOCK 1024

__global__ void histogramGPU ( unsigned int* data, unsigned int imageLength, unsigned int* histogram );
__global__ void saturate ( unsigned int *intHistogram, unsigned int imageLength );

int main ( int argc, char *argv[] ) {
	struct pb_Parameters *parameters = pb_ReadParameters ( &argc, argv );

	printf ( "Base implementation of histogramming.\n" );
	printf ( "Maintained by Nady Obeid <obeid1@ece.uiuc.edu>\n" );

	if ( !parameters ) {
		return -1;
	}

	if ( !parameters->inpFiles[0] ) {
		fputs ( "Input file expected\n", stderr );
		return -1;
	}

	int numberOfIterations = 0;
	if ( argc >= 2 ) {
		numberOfIterations = atoi ( argv[1] );
	} else {
		fputs ( "Expected at least one command line argument\n", stderr );
		return -1;
	}

	struct pb_TimerSet timers;
	pb_InitializeTimerSet ( &timers );

	char *inputStr = "Input";
	char *outputStr = "Output";

	pb_AddSubTimer ( &timers, inputStr, pb_TimerID_IO );
	pb_AddSubTimer ( &timers, outputStr, pb_TimerID_IO );

	pb_SwitchToSubTimer ( &timers, inputStr, pb_TimerID_IO );

	unsigned int img_width = 0;
	unsigned int img_height = 0;
	unsigned int histo_width = 0;
	unsigned int histo_height = 0;

	FILE *file = fopen ( parameters->inpFiles[0], "rb" );

	int result = 0;

	result += fread ( &img_width, sizeof(unsigned int), 1, file );
	result += fread ( &img_height, sizeof(unsigned int), 1, file );
	result += fread ( &histo_width, sizeof(unsigned int), 1, file );
	result += fread ( &histo_height, sizeof(unsigned int), 1, file );

	if ( result != 4 ) {
		fputs ( "Error reading input and output dimensions from file\n", stderr );
		return -1;
	}

	unsigned int *image = ( unsigned int* ) malloc ( img_width * img_height * sizeof ( unsigned int ) );
	unsigned int *histogram = ( unsigned int* ) calloc ( histo_width * histo_height, sizeof ( unsigned int ) );

	pb_SwitchToSubTimer ( &timers, "Input", pb_TimerID_IO );

	result = fread ( image, sizeof(unsigned int), img_width * img_height, file );

	fclose ( file );

	if ( result != img_width * img_height ) {
		fputs ( "Error reading input array from file\n", stderr );
		return -1;
	}

	unsigned int *deviceImage = NULL;
	unsigned int *intHistogram = NULL;

	cudaMalloc ( ( void** ) &deviceImage, img_width * img_height * sizeof ( unsigned int ) );
	cudaMalloc ( ( void** ) &intHistogram, histo_width * histo_height * sizeof ( unsigned int ) );

	cudaMemcpy ( deviceImage, image, img_width * img_height * sizeof ( unsigned int ), cudaMemcpyHostToDevice );

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	for ( int iteration = 0; iteration < numberOfIterations; ++iteration ) {
		unsigned int imageLength = img_width * img_height;
		unsigned int histogramLength = histo_width * histo_height;

		cudaMemset ( intHistogram, 0, histogramLength * sizeof ( unsigned int ) );

		dim3 imageBlockDimensions ( MAX_THREADS_PER_BLOCK, 1, 1 );
		int imageBlocks = ( imageLength + MAX_THREADS_PER_BLOCK - 1 ) / MAX_THREADS_PER_BLOCK;
		dim3 imageGridDimensions ( imageBlocks, 1, 1 );
		histogramGPU<<<imageBlockDimensions, imageGridDimensions>>> ( deviceImage, imageLength, intHistogram );

		dim3 histogramBlockDimensions ( MAX_THREADS_PER_BLOCK, 1, 1 );
		int histogramBlocks = ( histogramLength + MAX_THREADS_PER_BLOCK - 1 ) / MAX_THREADS_PER_BLOCK;
		dim3 histogramGridDimensions ( histogramBlocks, 1, 1 );
		saturate<<<histogramBlockDimensions, histogramGridDimensions>>> ( intHistogram, histogramLength );
	}

	//  pb_SwitchToTimer(&timers, pb_TimerID_IO);
	pb_SwitchToSubTimer ( &timers, outputStr, pb_TimerID_IO );

	cudaMemcpy ( histogram, intHistogram, histo_width * histo_height * sizeof ( unsigned int ), cudaMemcpyDeviceToHost );

	cudaFree ( deviceImage );
	cudaFree ( intHistogram );

	if ( parameters->outFile ) {
		dump_histo_img ( histogram, histo_height, histo_width, parameters->outFile );
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	free ( image );
	free ( histogram );

	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );

	printf ( "\n" );
	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );

	return 0;
}
