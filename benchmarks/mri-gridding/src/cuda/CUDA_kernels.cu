/***************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "UDTypes.h"

#define max(x, y) ( ( x < y ) ? y : x )
#define min(x, y) ( ( x > y ) ? y : x )

#define PI 3.14159265359

__device__ float device_kernel_value_CPU ( float v ) {
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

__device__ float device_kernel_value_LUT ( float v, float *LUT, int sizeLUT, float _1overCutoff2 ) {
    unsigned int k0;
    float v0;

    v *= (float)sizeLUT;
    k0 = (unsigned int)(v * _1overCutoff2);
    v0 = ((float)k0) / _1overCutoff2;
    return LUT[k0] + ((v - v0) * (LUT[k0 + 1] - LUT[k0]) / _1overCutoff2);
}

__global__ void gridding_Gold_kernel ( int n, Parameters *parameters, ReconstructionSample *samples, float *LUT, unsigned int sizeLUT, cmplx *gridData, float *sampleDensity ) {
    int i = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    if ( i >= n ) {
        return;
    }

    ReconstructionSample sample = samples[i];

    const unsigned int size_x = parameters->gridSize[0];
    const unsigned int size_y = parameters->gridSize[1];
    const unsigned int size_z = parameters->gridSize[2];

    const float cutoff         = ( ( float ) parameters->kernelWidth ) / 2.0;    
    const float cutoff2        = cutoff * cutoff;                         
    const float oneOverCutoff2 = 1 / cutoff2;                             

    const float beta = PI * sqrt ( 4 * parameters->kernelWidth * parameters->kernelWidth / ( parameters->oversample * parameters->oversample ) * ( parameters->oversample - .5 ) * ( parameters->oversample - .5 ) - .8 );

    const float sampleX = sample.kX;
    const float sampleY = sample.kY;
    const float sampleZ = sample.kZ;

    const unsigned int minX = max ( ( sampleX - cutoff ), 0.0 );
    const unsigned int maxX = min ( ( sampleX + cutoff ), size_x - 1.0 );
    const unsigned int minY = max ( ( sampleY - cutoff ), 0.0 );
    const unsigned int maxY = min ( ( sampleY + cutoff ), size_y - 1.0 );
    const unsigned int minZ = max ( ( sampleZ - cutoff ), 0.0 );
    const unsigned int maxZ = min ( ( sampleZ + cutoff ), size_z - 1.0 );

    const int useLUT = parameters->useLUT;

    for ( float z = minZ; z <= maxZ; ++z ) {
        for ( float y = minY; y <= maxY; ++y ) {
            for ( float x = minX; x <= maxX; ++x ) {

                /* value to evaluate kernel at */
                const float distanceZSquared = ( ( sampleZ - z ) * ( sampleZ - z ) );
                const float distanceYSquared = ( ( sampleY - y ) * ( sampleY - y ) );
                const float distanceXSquared = ( ( sampleX - x ) * ( sampleX - x ) );

                const float distance  = distanceZSquared + distanceYSquared + distanceXSquared;

                if ( distance < cutoff2 ) {
                    /* linear offset into 3-D matrix to get to zposition */
                    const unsigned int offsetZ = z * size_x * size_y;
                    const unsigned int offsetY = y * size_x;
                    const unsigned int offzetX = x;

                    /* linear index of (x,y,z) point */
                    const unsigned int index = offzetX + offsetY + offsetZ;

                    /* kernel weighting value */
                    float weight = 0;

                    if ( useLUT ) {
                        weight = device_kernel_value_LUT ( distance, LUT, sizeLUT, oneOverCutoff2 ) * sample.sdc;
                    } else {
                        weight = device_kernel_value_CPU ( beta * sqrt ( 1.0 - ( distance * oneOverCutoff2 ) ) ) * sample.sdc;
                    }

                    /* grid data */
                    // gridData[index].real += ( weight * sample.real );
                    // gridData[index].imag += ( weight * sample.imag );

                    float real = ( weight * sample.real );
                    float imag = ( weight * sample.imag );
                    atomicAdd ( &gridData[index].real, real );
                    atomicAdd ( &gridData[index].imag, imag );

                    /* estimate sample density */
                    // sampleDensity[index] += 1.0;
                    atomicAdd ( &sampleDensity[index], 1.0 );
                }
            }
        }
    }
}

#define MAX_NUMBER_OF_THREADS 1024

extern "C"
void gridding_Gold ( unsigned int n, Parameters parameters, ReconstructionSample *samples, float *LUT, unsigned int sizeLUT, cmplx *gridData, float *sampleDensity ) {
    const int gridNumberOfElements = parameters.gridSize[0] * parameters.gridSize[1] * parameters.gridSize[2];

    Parameters *deviceParameters = NULL;
    ReconstructionSample *deviceSamples = NULL;
    float *deviceLUT = NULL;
    cmplx *deviceGridData = NULL;
    float *deviceSampleDensity = NULL;

    cudaMalloc ( ( void** ) &deviceParameters, sizeof ( Parameters ) );
    cudaMalloc ( ( void** ) &deviceSamples, n * sizeof ( ReconstructionSample ) );
    cudaMalloc ( ( void** ) &deviceLUT, sizeLUT * sizeof ( float ) );
    cudaMalloc ( ( void** ) &deviceGridData, gridNumberOfElements * sizeof ( cmplx ) );
    cudaMalloc ( ( void** ) &deviceSampleDensity, gridNumberOfElements * sizeof ( float ) );

    cudaMemcpy ( deviceParameters, &parameters, sizeof ( Parameters ), cudaMemcpyHostToDevice );
    cudaMemcpy ( deviceSamples, samples, n * sizeof ( ReconstructionSample ), cudaMemcpyHostToDevice );
    cudaMemcpy ( deviceLUT, LUT, sizeLUT * sizeof ( float ), cudaMemcpyHostToDevice );

    cudaMemset ( deviceGridData, 0, gridNumberOfElements * sizeof ( cmplx ) );
    cudaMemset ( deviceSampleDensity, 0, gridNumberOfElements * sizeof ( float ) );

    // kernel call
    int numberOfBlocks = ( n + MAX_NUMBER_OF_THREADS - 1 ) / MAX_NUMBER_OF_THREADS;

    dim3 blockDimensions ( MAX_NUMBER_OF_THREADS, 1, 1 );
    dim3 gridDimensions ( numberOfBlocks, 1, 1 );

    gridding_Gold_kernel<<<gridDimensions, blockDimensions>>> ( n, deviceParameters, deviceSamples, deviceLUT, sizeLUT, deviceGridData, deviceSampleDensity );

    cudaMemcpy ( gridData, deviceGridData, gridNumberOfElements * sizeof ( cmplx ), cudaMemcpyDeviceToHost );
    cudaMemcpy ( sampleDensity, deviceSampleDensity, gridNumberOfElements * sizeof ( float ), cudaMemcpyDeviceToHost );

    cudaFree ( deviceParameters );
    cudaFree ( deviceSamples );
    cudaFree ( deviceLUT );
    cudaFree ( deviceGridData );
}