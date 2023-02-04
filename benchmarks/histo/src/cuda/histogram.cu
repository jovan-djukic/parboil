#define UINT8_MAX 255

__global__ void histogramGPU ( unsigned int* image, unsigned int imageLength, unsigned int* histogram ) {
    int tx = threadIdx.x; 
    int bx = blockIdx.x;

    int index = ( bx * blockDim.x ) + tx;

    if ( index < imageLength ) {
        unsigned int value = image[index];

        atomicAdd ( &histogram[value], 1 );
    }
}

__global__ void saturate ( unsigned int *intHistogram, unsigned int histogramLength ) {
    int tx = threadIdx.x; 
    int bx = blockIdx.x;
    int index = ( bx * blockDim.x ) + tx;

    if ( index < histogramLength ) {
        if ( intHistogram[index] > UINT8_MAX ) {
            intHistogram[index] = UINT8_MAX;
        }     
    }
}