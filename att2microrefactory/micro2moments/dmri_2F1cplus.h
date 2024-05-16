/*==========================================================
 * dmri_2F1mex.c
 *
 * Implements dmri_2F1mex as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "hypergeom2F1.h"
#include "mexToMathsTypes.h"
#include "threadHelper.h"

class ThArgs : public DMRIThreader
{
public:
    ElementType g;
    BufferType z;
    BufferType f;
};

THFCNRET dmri_2F1_process_fcn( void* );

int dmri_2F1cplus(double* output,const double gamma,const unsigned long N, double* x, const unsigned int maxNumberOfThreads);
