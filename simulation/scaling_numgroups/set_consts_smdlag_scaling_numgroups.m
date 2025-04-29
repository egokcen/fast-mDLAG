% set_consts_smdlag_scaling_numgroups.m
%
% Description: Define constants for smDLAG model fitting.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Setup directories/path

restoredefaultpath;
addpath(dataDir);
addpath(resultDir);
addpath(perfDir);
addpath(genpath('../../smDLAG'));

%% mDLAG fitting parameters

xDim_fit = xDim;   % In these experiments we'll assume the ground truth
nyqperiod = 80;    % Minimum sampling period (in ms) that we can get away
                   % with when choosing the number of inducing points.
                   % Empirically, for the squared exponential function, 
                   % when tau = 100 ms, nearly 100% of the power spectrum 
                   % is contained within 6 Hz (i.e., ~12 Hz Nyquist rate).
% We'll use as few inducing points as are needed to maintain the
% Nyquist rate for these data (easy to determine, since we know the
% ground truth timescales).
S = ceil(Ttotal * binWidth / nyqperiod);
randomSeed = 1;    % Seed the random number generator for reproducibility
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
learnTau = true;      % Toggle whether to learn timescale parameters
learnDelays = true;   % Toggle whether to learn delay parameters
learnInducingLocs = false;  % Toggle whether to learn inducing point locations
verbose = false;      % Print fitting progress
minVarFrac = 0;       % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = false;       % For speed-up, remove latents that become inactive in all groups
saveCcov = false;     % Set to false to save memory when saving final results
saveWcov = false;     % Set to false to save memory when saving final results
     