function options = create_micro2shodf_options(M, N, P,mask)
    % Initialize the options structure with default values
    options.L = 8;                      % Maximum order of SH
    options.lambda = 0.0005;             % Tikhonov regularization parameter
    options.tl = 1.0e-6;                % Lower threshold for signal intensity
    options.tu = 1-1.0e-6;        % Upper threshold for signal intensity
    options.ADC0 = 3.0e-3;              % Diffusivity of free water
    options.optimal = true;             % Use optimal (global) fitting approach
    options.mask = mask;       % Voxel mask to specify the region of interest
    
    % Microstructure model sanity check options
    options.chkmod = true;              % Perform sanity checks on microstructure model
    options.flperp = 0.001;             % Lower fraction limit for perpendicular diffusivity
    options.Flperp = 0.999;             % Upper fraction limit for perpendicular diffusivity
    
    % Additional processing options
    options.chunksz = 1000;             % Chunk size for processing to improve performance
    options.recrop = false;             % Recrop signal after free-water subtraction
    options.bth = 1;                    % Threshold for considering b-values as belonging to the same shell

    % You can add or modify any additional options as needed here
end
