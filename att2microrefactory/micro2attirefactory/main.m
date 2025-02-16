
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
sh = micro2shodf(atti, gi, bi, lpar, lperp, [], options);
