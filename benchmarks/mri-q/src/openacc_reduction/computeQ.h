#ifndef _COMPUTE_Q_H_
#define _COMPUTE_Q_H_

#define PI 3.1415926535897932384626433832795029f
#define PIx2 6.2831853071795864769252867665590058f

#define MIN(X, Y) ( (X) < (Y) ? (X) : (Y) )
#define K_ELEMS_PER_GRID 2048

struct kValues {
    float Kx;
    float Ky;
    float Kz;
    float PhiMag;
};

void ComputePhiMagCPU ( int numK, float *phiR, float *phiI, float *phiMag );
void ComputeQCPU ( int numK, int numX, struct kValues *kVals, float *x, float *y, float *z, float *Qr, float *Qi );
void createDataStructsCPU ( int numK, int numX, float **phiMag, float **Qr, float **Qi );

#endif
