% set_consts_smdlag_basic_demo.m
%
% Description: This script defines constants for the basic demonstration 
%              of mDLAG fitting via inducing variables (smDLAG).
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Setup directories/path

restoredefaultpath;
addpath(dataDir);
addpath(resultDir);
addpath(genpath('../../smDLAG'));

%% mDLAG fitting parameters

xDim_fit = 8;
S = 50;            % Number of inducing points (used 30 or 100)
covType = 'rbf';   % Type of covariance kernel
randomSeed = 0;    % Seed the random number generator for reproducibility
prior_val = 1e-12; % Set to very small value to set uninformative prior      
prior.d.beta = prior_val; 
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;

tol = 1e-8;           % Tolerance to determine fitting convergence
maxIters = 3e4;       % Maximum fitting iterations
freqLB = 1;           % Check for convergence of lower bound every freqLB iterations
freqParam = 1;        % Store intermediate delay and timescale estimates every freqParam iterations
learnTau = true;      % Toggle whether to learn timescale parameters
learnDelays = true;   % Toggle whether to learn delay parameters
learnInducingLocs = false;  % Toggle whether to learn inducing point locations
verbose = true;       % Print fitting progress
minVarFrac = 0;       % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = false;       % For speed-up, remove latents that become inactive in all groups
saveCcov = false;     % Set to false to save memory when saving final results
saveWcov = false;     % Set to false to save memory when saving final results
        