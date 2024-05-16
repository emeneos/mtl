function builddll(target, outputName)
    % Define entry points
    entryPoint1 = 'micro2moments';
    entryPoint2 = 'atti2micro';  
    entryPoint3 = 'micro2shodf';  

    % Configuration for generating C++ code as a DLL
    cfg = coder.config(target);
    cfg.TargetLang = 'C++';
    cfg.GenerateReport = true;
    cfg.LaunchReport = false;

    % Include directories and additional source files
    cfg.CustomInclude = {'D:\uvalladolid\matlab\labcode\att2microrefactory\micro2moments'};   
    cfg.CustomSource = {'threadHelper.cpp','dmri_2F1cplus.cpp','hyperGeom2F1.cpp','generateSHMatrix.cpp','sh2hot.cxx','sphericalHarmonics.cpp'};
    cfg.CustomSourceCode = ['#include "generateSHMatrix.h"'];

    % Load dataset
    load test_data.mat;
    
    % Determine the size of the 'atti' variable
    [M, N, P, G] = size(atti);

    % Process atti data to microstructure parameters
    optionsAtti2micro = create_atti2micro_options(mask, M, N, P, G);
    [lpar, lperp, f] = atti2micro(atti, gi, bi, optionsAtti2micro);

    % Process data for SH and DTI conversion
    optionsMicro2shodf = create_micro2shodf_options(M, N, P, mask);
    sh = micro2shodf(atti, gi, bi, lpar, lperp, f, optionsMicro2shodf);

    tens = shadc2dti(sh, 'mask', mask, 'unroll', false);
    u0 = dti2xyz(tens, 'mask', mask);
    options = create_micro2moments_options(mask, u0);
    
    % Define arguments for both entry-point functions
    args1 = {sh, lpar, lperp, [], options};
    args2 = {atti, gi, bi, optionsAtti2micro};  
    args3 = {atti, gi, bi, lpar, lperp, [], optionsMicro2shodf};

    % Check if a custom name is provided, use a default if not specified
    if nargin < 2 || isempty(outputName)
        outputName = 'mySharedLib';
    end

    % Generate standalone DLL library for the functions
    codegen('-config', cfg, '-o', outputName, entryPoint1, '-args', args1, entryPoint2, '-args', args2, entryPoint3, '-args', args3);
end
