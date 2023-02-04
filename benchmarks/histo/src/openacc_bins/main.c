/***************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ***************************************************************************/

#include <parboil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"

#define UINT8_MAX 255

#include <openacc.h>
#include <cuda.h>

/******************************************************************************
* Implementation: Reference
* Details:
* This implementations is a scalar, minimally optimized version. The only 
* optimization, which reduces the number of pointer chasing operations is the 
* use of a temporary pointer for each row.
******************************************************************************/

int main(int argc, char* argv[]) {
  struct pb_TimerSet timers;
  struct pb_Parameters *parameters;

  printf("OPENACC implementation of histogramming.\n");

  parameters = pb_ReadParameters(&argc, argv);
  if (!parameters)
    return -1;

  if(!parameters->inpFiles[0]){
    fputs("Input file expected\n", stderr);
    return -1;
  }

  int numberOfIterations;
  if (argc >= 2){
    numberOfIterations = atoi(argv[1]);
  } else {
    fputs("Expected at least one command line argument\n", stderr);
    return -1;
  }

  pb_InitializeTimerSet(&timers);
  
  char *inputStr = "Input";
  char *outputStr = "Output";
  
  pb_AddSubTimer(&timers, inputStr, pb_TimerID_IO);
  pb_AddSubTimer(&timers, outputStr, pb_TimerID_IO);
  
  pb_SwitchToSubTimer(&timers, inputStr, pb_TimerID_IO);  

  unsigned int imageWidth, imageHeight;
  unsigned int histogramWidth, histogramHeight;

  FILE* f = fopen(parameters->inpFiles[0],"rb");
  int result = 0;

  result += fread(&imageWidth,    sizeof(unsigned int), 1, f);
  result += fread(&imageHeight,   sizeof(unsigned int), 1, f);
  result += fread(&histogramWidth,  sizeof(unsigned int), 1, f);
  result += fread(&histogramHeight, sizeof(unsigned int), 1, f);

  if (result != 4){
    fputs("Error reading input and output dimensions from file\n", stderr);
    return -1;
  }

  unsigned int* image = (unsigned int*) malloc (imageWidth*imageHeight*sizeof(unsigned int));
  unsigned int* histogram = (unsigned int*) calloc (histogramWidth*histogramHeight, sizeof(unsigned int));
  
  pb_SwitchToSubTimer(&timers, "Input", pb_TimerID_IO);

  result = fread(image, sizeof(unsigned int), imageWidth*imageHeight, f);

  fclose(f);

  if (result != imageWidth*imageHeight){
    fputs("Error reading input array from file\n", stderr);
    return -1;
  }

  pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

  const int imageLength = imageWidth * imageHeight;
  const int histogramLength = histogramWidth * histogramHeight;

  const int vectorSize = 256;
  const int numberOfLocalHistograms = 2;
  const int localHistogramsLength = histogramLength * numberOfLocalHistograms;

  const int numberOfGangs = ( int ) ceil ( ( double ) imageLength / vectorSize );

  unsigned int *localHistograms = ( unsigned int* ) malloc ( sizeof ( unsigned int )  * localHistogramsLength );

  #pragma acc data copyin(image[0:imageLength]) copyout(histogram[0:histogramLength]) create(localHistograms[0:localHistogramsLength])
  {
    for (int iter = 0; iter < numberOfIterations; iter++){

      // #pragma acc parallel loop
      // for ( unsigned int i = 0; i < histogramLength; ++i ) {
      //   histogram[i] = 0;
      // }

      cudaMemset ( acc_deviceptr ( ( void* ) histogram ), 0, histogramLength * sizeof ( unsigned int ) );

      #pragma acc parallel loop  
      for ( unsigned int i = 0; i < localHistogramsLength; ++i ) {
        localHistograms[i] = 0;
      }

      #pragma acc parallel loop num_gangs(numberOfGangs) vector_length(vectorSize)
      for ( unsigned int i = 0; i < imageLength; ++i) {
        const unsigned int workersPerGang = ( unsigned int ) ceil ( ( double ) imageLength / numberOfGangs );
        const unsigned int gangId = i / workersPerGang;
        const unsigned int localHistogramNumber = gangId % numberOfLocalHistograms;

        const unsigned int value = image[i];

        #pragma acc atomic update
        ++localHistograms[histogramLength * localHistogramNumber + value];
      }

      #pragma acc parallel loop
      for ( unsigned int i = 0; i < localHistogramsLength; ++i ) {
        const unsigned int index = i % histogramLength;
        #pragma acc atomic update
        histogram[index] += localHistograms[i];
      }

      #pragma acc parallel loop
      for ( int i = 0; i < histogramLength; ++i ) {
        if ( histogram[i] > UINT8_MAX ) {
          histogram[i] = UINT8_MAX;
        }
      }
    }
  }

  // for (int iter = 0; iter < numberOfIterations; iter++){

  //   memset(histogram,0,histogramHeight*histogramWidth*sizeof(unsigned int));

  //   for (unsigned int i = 0; i < imageWidth*imageHeight; ++i) {
  //     const unsigned int value = image[i];
  //     // if (histogram[value] < UINT8_MAX) {
  //       ++histogram[value];
  //     // } 

  //   }
  // }


//  pb_SwitchToTimer(&timers, pb_TimerID_IO);
  pb_SwitchToSubTimer(&timers, outputStr, pb_TimerID_IO);

  if (parameters->outFile) {
    dump_histo_img(histogram, histogramHeight, histogramWidth, parameters->outFile);
  }

  pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

  free(image);
  free(histogram);

  pb_SwitchToTimer(&timers, pb_TimerID_NONE);

  printf("\n");
  pb_PrintTimerSet(&timers);
  pb_FreeParameters(parameters);

  return 0;
}
