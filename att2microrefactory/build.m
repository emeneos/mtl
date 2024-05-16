function build(target)
%build mex or build lib


    entryPoint = 'atti2micro';


    cfg = coder.config(target);

    cfg.TargetLang = 'C++';

    cfg.CustomSource = ['mexGenerateSHMatrix.cpp']; %using code from the function mexGenerateSHMatrix.cpp

    cfg.CustomInclude='D:\uvalladolid\DMRIMatlab\mexcode\mathsmex'; %adding path where the code from above is

    
    cfg.CustomSourceCode = ['#include "mexGenerateSHMatrix.h"']; %you need a header 
    
   
    cfg.GenerateReport = true;
    cfg.LaunchReport = false;

       
    load test_data.mat;
    T=toc; % This is always a large piece of data
    fprintf(1,'It took %f seconds to load data\n',T);
    
    tic;
    [M,N,P,G] = size(atti);
    options = create_atti2micro_options(mask,M,N,P,G); % you need a way to pass the data and other options to the line below
    
    
    %generate code
    codegen(entryPoint,'-args',{atti, gi,bi, options},'-config', cfg); % here the magic starts, you can modify {atti, gi,bi, options} and pass your own parameters here

end