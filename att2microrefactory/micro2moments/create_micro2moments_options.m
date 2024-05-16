function options = create_micro2moments_options(mask,u0)
    % Initialize the options structure with default values and specified parameters
                  % Double, either 1x1 or empty
    options.mask = mask;    % Boolean the same size as the image field
     if isempty(u0)
        options.u0 = [];             % Ensure u0 is an empty array if input u0 is empty
        options.type = 'rtpp';         % Variable length string
        options.nu = 0; 


    else
        options.u0 = u0;             % Double, either MxNxPx3 or empty
        options.type = 'ea';         % Variable length string
        options.nu = -1/2;  
    end
    
    % Microstructure model sanity check options
    options.chkmod = true;         % Always 1x1 boolean
    options.flperp = 0.001;        % Always 1x1 double
    options.Flperp = 0.999;        % Always 1x1 double
    options.ADC0 = 3.0e-3;         % Always 1x1 double
    
    % Additional processing options
    options.tau = 70.0e-3;         % Always 1x1 double
    options.chunksz = 256;         % Always 1x1 double
    options.clean = 0;             % Always 1x1 double
end