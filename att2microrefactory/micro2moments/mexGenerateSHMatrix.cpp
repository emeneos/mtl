/*==========================================================
 * mexGenerateSHMatrix.c
 *
 * Re-implements GenerateSHMatrix as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 * nlhs: Number of left-hand side arguments (outputs).
 * plhs: Array of pointers to mxArray pointers for left-hand side arguments (outputs). plhs[0] is the SH matrix
 * nrhs: Number of right-hand side arguments (inputs).
 * prhs: Array of pointers to mxArray pointers for right-hand side arguments (inputs).
 * in this case prhs[0] is L (Order of Spherical Harmonics) and prhs[1] is G (Gradient Directions or Angular Coordinates)
 *========================================================*/





#include "mex.h"
/*include <D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/>*/
#include "D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/sphericalHarmonics.h"
#include "D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/mexToMathsTypes.h"
#ifdef CODER 
#include "mex.h"
#include "D:\uvalladolid\matlab\labcode\att2microrefactory\micro2moments\sphericalHarmonics.cpp"
#include "D:\uvalladolid\matlab\labcode\att2microrefactory\micro2moments\sphericalHarmonics.h"
#include "D:\uvalladolid\matlab\labcode\att2microrefactory\micro2moments\mexToMathsTypes.h"
int generateSHMatrixx( double* plhs0, double* plhs1, const unsigned int L,  const double* Gi, const unsigned int G_ )
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
 /* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    /* check for proper number of arguments */
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:nrhs", "At least two inputs required.");
    }
    if (nlhs < 1) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:nlhs", "At least one output required.");
    }

   
    /* make sure the first input argument is scalar */
    if (!mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]) ||
        mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:notScalar", "Input L must be a scalar.");
    }

    /* make sure the first argument is an integer */
    double Ld = mxGetScalar(prhs[0]);
    unsigned int L = (unsigned int)Ld;
    double err = Ld - L;
    err = (err > 0 ? err : -err);
    if (err > 10 * mxGetEps()) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:notInteger", "Input L must be an integer.");
    }

    /* make sure the first argument is even */
    if (L != 2 * (L / 2)) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:notEven", "Input L must be even.");
    }

    /* make sure the second input argument is type double */
    if (!mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:notDouble", "Input G must be type double.");
    }

    /* make sure the second argument has proper size */
    size_t G = mxGetM(prhs[1]);
    size_t D = mxGetN(prhs[1]);
    if (D != 3) {
        mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:timesThree", "Input G must be Nx3.");
    }

    /* create a pointer to the real data in the input matrix  */
    mxDouble* Gi = mxGetDoubles(prhs[1]);

    double* gx = new double[G];
    double* gy = new double[G];
    double* gz = new double[G];
    for (unsigned int k = 0; k < G; ++k) {
        gx[k] = Gi[k];
        gy[k] = Gi[k + G];
        gz[k] = Gi[k + 2 * G];
    }

    /* if a second output is requested, compute the (unique) associated Legendre polynomials of degree L */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(G, L + 1, mxREAL);
        double* buffer = new double[L + 1];
        for (unsigned int k = 0; k < G; ++k) {
            shmaths::computeAssociatedLegendrePolynomialsL(gz[k], L, buffer);
            for (unsigned int l = 0; l <= L; ++l)
                mxGetDoubles(plhs[1])[l * G + k] = buffer[l];
        }
        delete[] buffer;
    }
    /* make sure the third input argument is provided */
    // if (nrhs < 3) {
    //     mexErrMsgIdAndTxt("MyToolbox:mexGenerateSHMatrix:nrhs", "Third input (output matrix) required.");
    // }

    double* theta = new double[G];
    double* phi = new double[G];

    shmaths::computeSphericalCoordsFromCartesian(gx, gy, gz, theta, phi, G);

    delete[] gx;
    delete[] gy;
    delete[] gz;

    /* Access the output matrix directly from prhs[2] */
    // mxDouble* outputData = mxGetDoubles(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(G, shmaths::getNumberOfEvenAssociatedLegendrePolynomials(L), mxREAL); 
    shmaths::computeSHMatrixSymmetric(G, theta, phi, L, mxGetDoubles(plhs[0]));
    /* Compute the SH matrix and store the result in the provided matrix */
    // shmaths::computeSHMatrixSymmetric(G, theta, phi, L, outputData);
    delete[] theta;
    delete[] phi;

}
#endif




