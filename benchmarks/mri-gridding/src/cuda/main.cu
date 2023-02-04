/***************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <parboil.h>

#include "UDTypes.h"
#include "CUDA_kernels.h"

#define PI 3.14159265

/************************************************************ 
 * This function reads the parameters from the file provided
 * as a comman line argument.
 ************************************************************/
void setParameters ( FILE *file, Parameters *p ) {
    fscanf ( file, "aquisition.numsamples=%d\n", &p->numSamples );
    fscanf ( file, "aquisition.kmax=%f %f %f\n", &p->kMax[0], &p->kMax[1], &p->kMax[2] );
    fscanf ( file, "aquisition.matrixSize=%d %d %d\n", &p->aquisitionMatrixSize[0], &p->aquisitionMatrixSize[1], &p->aquisitionMatrixSize[2] );
    fscanf ( file, "reconstruction.matrixSize=%d %d %d\n", &p->reconstructionMatrixSize[0], &p->reconstructionMatrixSize[1], &p->reconstructionMatrixSize[2] );
    fscanf ( file, "gridding.matrixSize=%d %d %d\n", &p->gridSize[0], &p->gridSize[1], &p->gridSize[2] );
    fscanf ( file, "gridding.oversampling=%f\n", &p->oversample );
    fscanf ( file, "kernel.width=%f\n", &p->kernelWidth );
    fscanf ( file, "kernel.useLUT=%d\n", &p->useLUT );

    printf ( "  Number of samples = %d\n", p->numSamples );
    printf ( "  Grid Size = %dx%dx%d\n", p->gridSize[0], p->gridSize[1], p->gridSize[2] );
    printf ( "  Input Matrix Size = %dx%dx%d\n", p->aquisitionMatrixSize[0], p->aquisitionMatrixSize[1], p->aquisitionMatrixSize[2] );
    printf ( "  Recon Matrix Size = %dx%dx%d\n", p->reconstructionMatrixSize[0], p->reconstructionMatrixSize[1], p->reconstructionMatrixSize[2] );
    printf ( "  Kernel Width = %f\n", p->kernelWidth );
    printf ( "  KMax = %.2f %.2f %.2f\n", p->kMax[0], p->kMax[1], p->kMax[2] );
    printf ( "  Oversampling = %f\n", p->oversample );
    printf ( "  GPU Binsize = %d\n", p->binsize );
    printf ( "  Use LUT = %s\n", (p->useLUT ) ? "Yes" : "No" );
}

/************************************************************ 
 * This function reads the sample point data from the kspace
 * and klocation files (and sdc file if provided) into the
 * sample array.
 * Returns the number of samples read successfully.
 ************************************************************/
unsigned int readSampleData ( Parameters params, FILE *uksdata_f, ReconstructionSample *samples ) {
    int count = 0;
    for ( unsigned int i = 0; i < params.numSamples; ++i ) {
        if ( feof ( uksdata_f ) ) {
            break;
        }

        fread ( ( (void*) &samples[i] ), sizeof ( ReconstructionSample ), 1, uksdata_f );
        count++;
    }

    float kScale[3];
    kScale[0] = ( (float) params.aquisitionMatrixSize[0] ) / ( ( (float) params.reconstructionMatrixSize[0] ) * ( (float) params.kMax[0] ) );
    kScale[1] = ( (float) params.aquisitionMatrixSize[1] ) / ( ( (float) params.reconstructionMatrixSize[1] ) * ( (float) params.kMax[1] ) );
    kScale[2] = ( (float) params.aquisitionMatrixSize[2] ) / ( ( (float) params.reconstructionMatrixSize[2] ) * ( (float) params.kMax[2] ) );

    int size_x = params.gridSize[0];
    int size_y = params.gridSize[1];
    int size_z = params.gridSize[2];

    float ax = ( kScale[0] * ( size_x - 1 ) ) / 2.0;
    float bx = (float) ( size_x - 1 ) / 2.0;

    float ay = ( kScale[1] * ( size_y - 1 ) ) / 2.0;
    float by = (float) ( size_y - 1 ) / 2.0;

    float az = ( kScale[2] * ( size_z - 1 ) ) / 2.0;
    float bz = (float) ( size_z - 1 ) / 2.0;

    for ( int n = 0; n < count; n++ ) {
        samples[n].kX = floor ( ( samples[n].kX * ax ) + bx );
        samples[n].kY = floor ( ( samples[n].kY * ay ) + by );
        samples[n].kZ = floor ( ( samples[n].kZ * az ) + bz );
    }

    return count;
}

float kernel_value_CPU ( float v ) {
    const float z = v * v;

    // polynomials taken from http://ccrma.stanford.edu/CCRMA/Courses/422/projects/kbd/kbdwindow.cpp
    float num = (  z * ( z * ( z * ( z * ( z * ( z * ( z * ( z * ( z * ( z * ( z * ( z * ( z *
                (  z * 0.210580722890567e-22f + 0.380715242345326e-19f ) +
                0.479440257548300e-16f ) + 0.435125971262668e-13f ) + 0.300931127112960e-10f ) +
                0.160224679395361e-7f ) + 0.654858370096785e-5f ) + 0.202591084143397e-2f ) +
                0.463076284721000e0f ) + 0.754337328948189e2f ) + 0.830792541809429e4f ) +
                0.571661130563785e6f ) + 0.216415572361227e8f ) + 0.356644482244025e9f ) +
                0.144048298227235e10f );

    float den = ( z * ( z * ( z - 0.307646912682801e4f ) + 0.347626332405882e7f ) - 0.144048298227235e10f );

    float rValue = -num / den;

    return rValue;
}

void calculateLUT ( float beta, float width, float **LUT, unsigned int *sizeLUT ) {
    const float cutoff2 = ( width * width ) / 4.0;

    if ( width > 0 ) {
        // compute size of LUT based on kernel width
        unsigned int size = (unsigned int) ( 10000 * width );

        // allocate memory
        (*LUT) = (float*) malloc ( size * sizeof ( float ) );

        for ( unsigned int k = 0; k < size; ++k ) {
            // compute value to evaluate kernel at
            // v in the range 0:(_width/2)^2
            float v = ( ((float) k) / ((float) size) ) * cutoff2;

            // compute kernel value and store
            (*LUT)[k] = kernel_value_CPU ( beta * sqrt ( 1.0 - ( v / cutoff2 ) ) );
        }

        (*sizeLUT) = size;
    }
}


int main ( int argc, char *argv[] ) {
    struct pb_Parameters *prms;
    struct pb_TimerSet timers;

    prms = pb_ReadParameters ( &argc, argv );
    pb_InitializeTimerSet ( &timers );

    pb_SwitchToTimer ( &timers, pb_TimerID_NONE );

    char uksdata[250];
    Parameters params;

    FILE *uksfile_f = NULL;
    FILE *uksdata_f = NULL;

    strcpy ( uksdata, prms->inpFiles[0] );
    strcat ( uksdata, ".data" );

    uksfile_f = fopen ( prms->inpFiles[0], "r" );
    if ( uksfile_f == NULL ) {
        printf ( "ERROR: Could not open %s\n", prms->inpFiles[0] );
        exit ( 1 );
    }

    printf ( "\nReading parameters\n" );

    if ( argc >= 2 ) {
        params.binsize = atoi ( argv[1] );
    } else { //default binsize value;
        params.binsize = 128;
    }

    setParameters ( uksfile_f, &params );

    pb_SwitchToTimer ( &timers, pb_TimerID_IO );

    ReconstructionSample *samples = (ReconstructionSample*) malloc ( params.numSamples * sizeof ( ReconstructionSample ) ); //Input Data
    float *LUT;                                                                                                       //use look-up table for faster execution on CPU (intermediate data)
    unsigned int sizeLUT;                                                                                             //set in the function calculateLUT (intermediate data)

    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];

    cmplx *gridData = (cmplx*) calloc ( gridNumElems, sizeof ( cmplx ) );      //Output Data
    float *sampleDensity = (float*) calloc ( gridNumElems, sizeof ( float ) ); //Output Data

    if ( samples == NULL ) {
        printf ( "ERROR: Unable to allocate memory for input data\n" );
        exit ( 1 );
    }

    if ( sampleDensity == NULL || gridData == NULL ) {
        printf ( "ERROR: Unable to allocate memory for output data\n" );
        exit ( 1 );
    }

    uksdata_f = fopen ( uksdata, "rb" );

    if ( uksdata_f == NULL ) {
        printf ( "ERROR: Could not open data file\n" );
        exit ( 1 );
    }

    printf ( "Reading input data from files\n" );

    unsigned int n = readSampleData ( params, uksdata_f, samples );
    fclose ( uksdata_f );

    if ( params.useLUT ) {
        printf ( "Generating Look-Up Table\n" );
        float beta = PI * sqrt ( 4 * params.kernelWidth * params.kernelWidth / ( params.oversample * params.oversample ) * ( params.oversample - .5 ) * ( params.oversample - .5 ) - .8 );
        calculateLUT ( beta, params.kernelWidth, &LUT, &sizeLUT );
    }

    pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

    gridding_Gold ( n, params, samples, LUT, sizeLUT, gridData, sampleDensity );

    pb_SwitchToTimer ( &timers, pb_TimerID_IO );

    int passed = 1;

    FILE *outfile;
    if ( !( outfile = fopen ( prms->outFile, "w" ) ) ) {
        printf ( "Cannot open output file!\n" );
    } else {
        fwrite ( &passed, sizeof ( int ), 1, outfile );
        fclose ( outfile );
    }

    pb_SwitchToTimer ( &timers, pb_TimerID_NONE );

    if ( params.useLUT ) {
        free ( LUT );
    }
    
    free ( samples );
    free ( gridData );
    free ( sampleDensity );

    printf ( "\n" );
    pb_PrintTimerSet ( &timers );
    pb_FreeParameters ( prms );

    return 0;
}
