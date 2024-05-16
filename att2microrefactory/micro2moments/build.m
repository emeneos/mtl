function build(target)
%build mex or build lib


    entryPoint = 'micro2moments';


    cfg = coder.config(target);

    cfg.TargetLang = 'C++';

    cfg.CustomSource = ['mexGenerateSHMatrix.cpp']; %using code from the function mexGenerateSHMatrix.cpp

    cfg.CustomInclude='D:\uvalladolid\DMRIMatlab\mexcode\mathsmex'; %adding path where the code from above is

    
    cfg.CustomSourceCode = ['#include "mexGenerateSHMatrix.h"']; %you need a header 
    
   
    cfg.GenerateReport = true;
    cfg.LaunchReport = false;

       
    % Load the dataset from a MAT file
    load test_data.mat;
    
    % Determine the size of the 'atti' variable to pass as parameters
    [M, N, P, G] = size(atti);
    
    % Create options for the atti2micro conversion
    options = create_atti2micro_options(mask, M, N, P, G);
    
    % Convert atti data to microstructure parameters using the simplest use case
    [lpar, lperp, f] = atti2micro(atti, gi, bi, options);
    

    % Create options structure
    options = create_micro2shodf_options(M,N,P,mask);
    
    % Call micro2shodf with the options structure
    sh = micro2shodf(atti, gi, bi, lpar, lperp, f, options);
    
    tic;
    
    tens = shadc2dti(sh,'mask',mask,'unroll',false ); % Compute a tensor field from the SH volume
    u0 = dti2xyz( tens, 'mask', mask ); % Maximum diffusion direction
    % Assuming 'mask' is already defined in your script
    options = create_micro2moments_options(mask,u0);
    
    
    

    
    %generate code
    codegen(entryPoint,'-args',{ sh, lpar, lperp, [], options},'-config', cfg); % here the magic starts, you can modify {atti, gi,bi, options} and pass your own parameters here

end