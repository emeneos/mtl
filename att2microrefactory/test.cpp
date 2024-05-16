


#include "mex.h"
#include "D:\uvalladolid\DMRIMatlab\mexcode\mathsmex\sphericalHarmonics.cpp"
#include "D:\uvalladolid\DMRIMatlab\mexcode\mathsmex\sphericalHarmonics.h"
#include "D:\uvalladolid\DMRIMatlab\mexcode\mathsmex\mexToMathsTypes.h"

#ifdef CODER 
int test( double* plhs0, double* plhs1, const unsigned int L,  const double* Gi, const unsigned int G_ )
{
    /* make sure the first argument is even */
    if (L != 2 * (L / 2)) {
        return -1;
    }

    size_t G = (size_t)G_;
    size_t D = 3;

    double* gx = new double[G];
    double* gy = new double[G];
    double* gz = new double[G];
    for (unsigned int k = 0; k < G; ++k) {
        gx[k] = Gi[k];
        gy[k] = Gi[k + G];
        gz[k] = Gi[k + 2 * G];
    }

    if(plhs1!=NULL){
        double* buffer = new double[L + 1];
        for (unsigned int k = 0; k < G; ++k) {
            shmaths::computeAssociatedLegendrePolynomialsL(gz[k], L, buffer);
            for (unsigned int l = 0; l <= L; ++l)
                plhs1[l * G + k] = buffer[l];
        }
        delete[] buffer;
    }

    double* theta = new double[G];
    double* phi = new double[G];

    shmaths::computeSphericalCoordsFromCartesian(gx, gy, gz, theta, phi, G);

    delete[] gx;
    delete[] gy;
    delete[] gz;

    shmaths::computeSHMatrixSymmetric(G, theta, phi, L, plhs0 );
    /* Compute the SH matrix and store the result in the provided matrix */
    // shmaths::computeSHMatrixSymmetric(G, theta, phi, L, outputData);
    delete[] theta;
    delete[] phi;

    return 0;

}

#else
int test( double* plhs0, double* plhs1, const unsigned int L,  const double* Gi, const unsigned int G_ )
{return 44;
}
#endif