
% Load the dataset from a MAT file
load test_data.mat;

% Measure and print the time taken to load the data
T = toc;
fprintf(1, 'It took %f seconds to load data\n', T);

% Start measuring time for processing
tic;

% Determine the size of the 'atti' variable to pass as parameters
[M, N, P, G] = size(atti);

% Create options for the atti2micro conversion
options = create_atti2micro_options(mask, M, N, P, G);

% Convert atti data to microstructure parameters using the simplest use case
[lpar, lperp, f] = atti2micro(atti, gi, bi, options);

% Measure and print the time taken for processing
T = toc;
fprintf(1, 'It took %f seconds to complete\n', T);


% Create options structure
options = create_micro2shodf_options(M,N,P,mask);

% Call micro2shodf with the options structure
sh = micro2shodf(atti, gi, bi, lpar, lperp, f, options);

tic;

tens = shadc2dti(sh,'mask',mask,'unroll',false ); % Compute a tensor field from the SH volume
u0 = dti2xyz( tens, 'mask', mask ); % Maximum diffusion direction
% Assuming 'mask' is already defined in your script
options = create_micro2moments_options(mask,u0);





% Then call micro2moments with these options
UEm2 = micro2moments( sh, lpar, lperp, [], options ); % Orders >-3 are allowed
T = toc; % For integer nu, the computation is quite fast
fprintf(1,'It took %f seconds to compute the moment\n',T);