% smdlag_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the smDLAG
%              model. Ground truth data and parameters are saved to a file.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     28 Nov 2022 -- Initial full revision.

%% Define smDLAG ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 100;                          % Number of sequences (trials)
T = 50;                           % Number of samples per sequence
S = 25;                           % Number of inducing points
binWidth = 20;                    % Sample period of ground truth data
yDims = [10 10 10];               % Dimensionalities of each observed group
yDim = sum(yDims);                % Joint dimensionality of all groups
M = length(yDims);                % Total number of groups
xDim = 7;                         % Across-group latent dimensionality
snr = 0.2*ones(1,M);              % Signal-to-noise ratios of each group

% a_alpha (M x xDim array) is where the inter-group interaction structure 
% is defined. Row i corresponds to group i. Column j corresponds to latent 
% j. A value of Inf indicates that a latent is NOT present in a group. The
% corresponding loadings will be 0 for that group.
MAG = 20; % Control the variance of alpha parameters (larger = less var.)
hyperparams.a_alpha = MAG.*[1 1   1   Inf 1   Inf Inf;
                            1 1   Inf 1   Inf 1   Inf;
                            1 Inf 1   1   Inf Inf 1  ];
% The remaining hyperparameters are not very important, and can be left
% alone.
hyperparams.b_alpha = MAG.*ones(M,xDim);
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
hyperparams.beta = 1;

% Gaussian process (GP) parameters
tau = [20 30 50 80 100 130 150];  % GP timescales
eps = 1e-3.*ones(1,xDim);         % GP noise variances
D = [0  0  0  0  0  0  0;         % Latent delay matrix
     15 30 0  0  0  0  0;
     30 0  25 40 0  0  0];
Z = repmat({linspace(1,T,S).*binWidth},xDim,1); % Inducing point locations

%% Randomly generate data from a mDLAG model

[seqTrue, paramsTrue] = simdata_smdlag(N, T, binWidth, yDims, xDim, ...
    hyperparams, snr, tau, eps, D, Z);

%% Visualize the ground truth

% GP parameters
sigDimsTrue = (1./hyperparams.a_alpha) > 0;
units = 'ms';
gp_params = plotGPparams_smdlag(paramsTrue,binWidth, ...
    'sigDims',sigDimsTrue,'units',units);

% Latent timecourses
% Relative shared variance explained by each dimension
alpha_inv = 1./paramsTrue.alphas;
alpha_inv_rel = alpha_inv ./ sum(alpha_inv,2);
[seqTrue, sortParams] = scaleByVarExp(seqTrue, paramsTrue, alpha_inv_rel, ...
                                      'sortDims', false, ...
                                      'sortGroup', 1, ...
                                      'numDim', 7, ...
                                      'indatafield', 'xsm', ...
                                      'outdatafield', 'xve');
xspec = 'xsm';
plotDimsVsTime_smdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', true, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {}, ...
                     'plotW', true);
                 
%% Save generated ground truth data and parameters

save('demo/data/smdlag_demo_data.mat', ...
    'seqTrue', 'paramsTrue', 'binWidth', 'sigDimsTrue');