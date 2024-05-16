% Start measuring time for data loading
tic;

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

% Define slices to be visualized
sl = [5, 7, 9, 11, 13, 15];

% Arrange the slices for 'f' parameter visualization
% Here, 'f' represents parallel diffusivity
F = [ f(:,:,sl(1))', f(:,:,sl(2))'; f(:,:,sl(3))', f(:,:,sl(4))'; f(:,:,sl(5))', f(:,:,sl(6))'];

% Arrange the slices for 'lpar' parameter visualization
% 'lpar' represents longitudinal (parallel) diffusivity
LPA = [ lpar(:,:,sl(1))', lpar(:,:,sl(2))'; lpar(:,:,sl(3))', lpar(:,:,sl(4))'; lpar(:,:,sl(5))', lpar(:,:,sl(6))'];

% Arrange the slices for 'lperp' parameter visualization
% 'lperp' represents perpendicular diffusivity
LPP = [ lperp(:,:,sl(1))', lperp(:,:,sl(2))'; lperp(:,:,sl(3))', lperp(:,:,sl(4))'; lperp(:,:,sl(5))', lperp(:,:,sl(6))'];

% Set the colormap based on the 'mapType' variable
mapType = 'gray'; % Default to grayscale
switch(mapType)
    case 'high'
        MAP = psychedelia(512); % Psychedelic colormap
    case 'gray'
        MAP = gray(512); % Grayscale colormap
    case 'default'
        MAP = parula(512); % Default MATLAB colormap
end

% Close the previous figure and create a new figure for visualization
close(figure(1));
hf = figure(1);
set(hf, 'Name', 'Micro-Structure Model', 'Position', [10, 10, 700, 340]);

% Visualization setup
R = 1; C = 3; % Rows and columns for subplot arrangement

% Display the 'f' parameter
r = 1; c = 1;
subplot('Position', [(c-1)/C, (r-1)/R, 1/C, 1/R]);
imagesc(F, [0, 1]);
colormap(MAP);
colorbar;
title('f (non-free water)');
axis('equal'); axis('off'); axis('tight');

% Display the 'lpar' parameter
r = 1; c = 2;
subplot('Position', [(c-1)/C, (r-1)/R, 1/C, 1/R]);
imagesc(LPA, [1.5e-3, 3.0e-3]);
colormap(MAP);
colorbar;
title('\lambda_{||} (mm^2/s)');
axis('equal'); axis('off'); axis('tight');

% Display the 'lperp' parameter
r = 1; c = 3;
subplot('Position', [(c-1)/C, (r-1)/R, 1/C, 1/R]);
imagesc(LPP, [0, 0.3e-3]);
colormap(MAP);
colorbar;
title('\lambda_{\perp} (mm^2/s)');
axis('equal'); axis('off'); axis('tight');
