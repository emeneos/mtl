function build(target)
%build mex or build lib


    entryPoint = 'micro2shodf';


    % Configuration for generating C++ code as a DLL
    cfg = coder.config('dll');
    cfg.TargetLang = 'C++';
    cfg.GenerateReport = true;
    cfg.LaunchReport = false;

    
    cfg.CustomInclude = ['/media/sf_att2microrefactory/micro2moments'];   
    % Specify additional source files and include directories
    cfg.CustomSource = ['threadHelper.cpp','dmri_2F1cplus.cpp','hyperGeom2F1.cpp','generateSHMatrix.cpp','sh2hot.cxx','sphericalHarmonics.cpp'];
    cfg.CustomSourceCode = ['#include "generateSHMatrix.h"'];

       
    load test_data.mat;
    T=toc; % This is always a large piece of data
    fprintf(1,'It took %f seconds to load data\n',T);
    
    tic;
    [M,N,P,G] = size(atti);
    options = create_atti2micro_options(mask,M,N,P,G); % you need a way to pass the data and other options to the line below
    % Convert atti data to microstructure parameters using the simplest use case
    [lpar, lperp, f] = atti2micro(atti, gi, bi, options);
    
    % Measure and print the time taken for processing
    T = toc;
    fprintf(1, 'It took %f seconds to complete\n', T);
    
    
    % Create options structure
    options = create_micro2shodf_options(M,N,P,mask);
    
    %generate code
    codegen(entryPoint,'-args',{atti, gi, bi, lpar, lperp, [], options},'-config', cfg); % here the magic starts, you can modify {atti, gi,bi, options} and pass your own parameters here

end