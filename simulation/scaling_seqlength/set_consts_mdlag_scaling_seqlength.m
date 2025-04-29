% set_consts_mdlag_scaling_seqlength.m
%
% Description: Define constants for mDLAG model fitting.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Setup directories/path

restoredefaultpath;
addpath(dataDir);
addpath(resultDir);
addpath(perfDir);
addpath(genpath('../../mDLAG'));

%% mDLAG fitting parameters

xDim_fit = xDim;   % In these experiments we'll assume the ground truth
randomSeed = 0;    % Seed the random number generator for reproducibility
prior_val = 1e-12; % Set to very small value to set uninformative prior      
prior.d.beta = prior_val; 
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;

tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 20000;     % Maximum fitting iterations
freqLB = 1;           % Check for convergence of lower bound every freqLB iterations
freqParam = 1;        % Store intermediate delay and timescale estimates every freqParam iterations
learnDelays = true;   % Toggle whether to learn delay parameters
verbose = false;      % Print fitting progress
minVarFrac = 0;       % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = false;       % For speed-up, remove latents that become inactive in all groups
saveXcov = false;     % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results
