#ifndef MEXGENERATESHMATRIX_H
#define MEXGENERATESHMATRIX_H

#include "mex.h"
/*#include <D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/>*/
#include "D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/sphericalHarmonics.h"
#include "D:/uvalladolid/DMRIMatlab/mexcode/mathsmex/mexToMathsTypes.h"

/* The gateway function */
void mexGenerateSHMatrix(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[]);
int generateSHMatrix( double* plhs0, double* plhs1, const unsigned int L,  const double* Gi, const unsigned int G_ );
#endif /* MEXGENERATESHMATRIX_H */