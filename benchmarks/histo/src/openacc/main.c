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

  int numIterations;
  if (argc >= 2){
    numIterations = atoi(argv[1]);
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

  unsigned int img_width, img_height;
  unsigned int histo_width, histo_height;

  FILE* f = fopen(parameters->inpFiles[0],"rb");
  int result = 0;

  result += fread(&img_width,    sizeof(unsigned int), 1, f);
  result += fread(&img_height,   sizeof(unsigned int), 1, f);
  result += fread(&histo_width,  sizeof(unsigned int), 1, f);
  result += fread(&histo_height, sizeof(unsigned int), 1, f);

  if (result != 4){
    fputs("Error reading input and output dimensions from file\n", stderr);
    return -1;
  }

  unsigned int* img = (unsigned int*) malloc (img_width*img_height*sizeof(unsigned int));
  unsigned int* histo = (unsigned int*) calloc (histo_width*histo_height, sizeof(unsigned int));
  
  pb_SwitchToSubTimer(&timers, "Input", pb_TimerID_IO);

  result = fread(img, sizeof(unsigned int), img_width*img_height, f);

  fclose(f);

  if (result != img_width*img_height){
    fputs("Error reading input array from file\n", stderr);
    return -1;
  }

  pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

  const int imageLength = img_width * img_height;
  const int histoLength = histo_width * histo_height;

  #pragma acc data copyin(img[0:imageLength]) copyout(histo[0:histoLength]) 
  {
    for (int iter = 0; iter < numIterations; iter++){

      // #pragma acc parallel loop
      // for ( unsigned int i = 0; i < histo_width * histo_height; ++i ) {
      //   histo[i] = 0;
      // }

      cudaMemset ( acc_deviceptr ( ( void* ) histo ), 0, histo_width * histo_height * sizeof ( unsigned int ) );
      #pragma acc parallel loop
      for ( unsigned int i = 0; i < img_width*img_height; ++i) {
        const unsigned int value = img[i];
        #pragma acc atomic update
        ++histo[value];
      }

      #pragma acc parallel loop
      for (unsigned int i = 0; i < histo_width * histo_height; ++i ) {
        if ( histo[i] > UINT8_MAX ) {
          histo[i] = UINT8_MAX;
        }
      }
    }
  }

  // for (int iter = 0; iter < numIterations; iter++){

  //   memset(histo,0,histo_height*histo_width*sizeof(unsigned int));

  //   for (unsigned int i = 0; i < img_width*img_height; ++i) {
  //     const unsigned int value = img[i];
  //     // if (histo[value] < UINT8_MAX) {
  //       ++histo[value];
  //     // } 

  //   }
  // }


//  pb_SwitchToTimer(&timers, pb_TimerID_IO);
  pb_SwitchToSubTimer(&timers, outputStr, pb_TimerID_IO);

  if (parameters->outFile) {
    dump_histo_img(histo, histo_height, histo_width, parameters->outFile);
  }

  pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

  free(img);
  free(histo);

  pb_SwitchToTimer(&timers, pb_TimerID_NONE);

  printf("\n");
  pb_PrintTimerSet(&timers);
  pb_FreeParameters(parameters);

  return 0;
}
