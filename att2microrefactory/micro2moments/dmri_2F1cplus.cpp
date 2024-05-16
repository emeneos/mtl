/*==========================================================
 * dmri_2F1cplus.cpp
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
int dmri_2F1cplus(double* output,const double gamma,const unsigned long N,  double* x, const unsigned int maxNumberOfThreads)
{
    ElementType g = gamma;

    unsigned int maxthreads = 1000000;
    if( maxNumberOfThreads>0 ){
        maxthreads = maxNumberOfThreads;
    }
    maxthreads = get_number_of_threads( maxthreads );

    //=======================================================================================
    unsigned int chunksz = (N/maxthreads)/20;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.g = g;        //
    threader.z = x; //
    threader.f = output; //
    //=======================================================================================
    threader.threadedProcess( maxthreads, dmri_2F1_process_fcn );
    //=======================================================================================
    return 1; // how to check if this was okay
}

THFCNRET dmri_2F1_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    ElementType g = args->g;
    
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            ElementType z = args->z[i];
            ElementType val;
            int result = hypegeo::hyperGeom2F1( g, z, val );
            args->f[i] = val;
        }
    }
    while( start < args->getN() );
    
    return (THFCNRET)NULL;
}
