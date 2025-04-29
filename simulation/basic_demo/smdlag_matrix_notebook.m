% smdlag_matrix_notebook.m
%
% Description: Visualize the various smDLAG GP covariance matrices.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Setup directories/path

restoredefaultpath;
addpath(genpath('../../smDLAG'));

%% Define smDLAG ground truth model parameters

% Set smDLAG model parameters. We only need to generate one trial for
% demonstration.
randomSeed = 0;                   % Seed for the random number generator
N = 1;                            % Total number of trials
T = 25;                           % Number of samples per sequence
S = 10;                           % Number of inducing points
binWidth = 20;                    % Sample period of ground truth data
yDims = [1 1];                    % Dimensionalities of each observed group
yDim = sum(yDims);                % Joint dimensionality of all groups
numGroups = length(yDims);        % Total number of groups
xDim = 1;                         % Across-group latent dimensionality
snr = ones(1,numGroups);          % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel


% Gaussian process (GP) parameters
tau = 100;                        % GP timescales
eps = 1e-5.*ones(1,xDim);         % GP noise variances
D = [0; 80];                      % GP time delays
Z = repmat({linspace(1,T,S).*binWidth},xDim,1); % Inducing point locations

%% Generate the various smDLAG covariance matrices

params = generate_params_smdlag(yDims, xDim, binWidth, ...
                                hyperparams, snr, tau, eps, D, Z);

Kx = construct_Kx_smdlag(params, T); Kx = Kx{1};
Kxw = construct_Kxw_smdlag(params,T); Kxw = Kxw{1};
Kw = construct_Kw_smdlag(params); Kw = Kw{1};

% Plot all covariance matrices
figure;
% Kx
subplot(1,3,1);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(Kx));
axis square;
axis off;
title('K^x');

% Kxw
subplot(1,3,2);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(Kxw));
axis off;
title('K^{xw}');

% Kw
subplot(1,3,3);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(Kw));
axis square;
axis off;
title('K^w');
