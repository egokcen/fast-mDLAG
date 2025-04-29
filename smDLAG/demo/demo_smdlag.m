% ===========================================================
% smDLAG DEMO: Run this script from the main smDLAG directory
% ===========================================================
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% ================
% 0) Load demo data
% =================

% Synthetic data generated from a smDLAG model.
% See simulation/smdlag_groundtruth_data.m for a demonstration of how to
% generate simulated data from a smDLAG model.
dat_file = 'demo/data/smdlag_demo_data';
fprintf('Reading from %s \n',dat_file);
load(dat_file);
numGroups = length(paramsTrue.yDims);
yDims = paramsTrue.yDims;
units = 'ms';

% Names for each group in the demo data file
groupNames = {'A', 'B', 'C'};

%% ======================================
% 1a) Initialize smDLAG model parameters
% ======================================

xDim_fit = 10;           % Number of latent variables to be fitted
S = 50;                  % Number of inducing point locations
randomSeed = 0;          % Seed the random number generator for reproducibility
prior_val = 1e-12;       % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;
saveCcov = false;        % Set to false to save memory when saving initParams
initParams = init_smdlag(seqTrue,yDims,xDim_fit,binWidth,...
                         'S', S, ...
                         'prior', prior, ...
                         'randomSeed', randomSeed, ...
                         'saveCcov', saveCcov);

%% ==========================
% 1b) Set fitting parameters
% ==========================

% Let's explicitly define relevant optional arguments, for
% the sake of demonstration:
tol = 1e-8;           % Tolerance to determine EM convergence
maxIters = 5000;      % Maximum EM iterations
freqLB = 1;           % Check for convergence of lower bound every freqLB EM iterations
freqParam = 1;        % Store intermediate delay and timescale estimates every freqParam EM iterations
learnTau = true;      % Toggle whether to learn timescale parameters
learnDelays = true;   % Toggle whether to learn delay parameters
learnInducingLocs = false; % Toggle whether to learn inducing point locations
verbose = true;       % Print fitting progress
minVarFrac = 0.001;   % Private noise variances will not be allowed to go below this value
maxDelayFrac  = 0.5;  % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveCcov = false;     % Set to false to save memory when saving final results
saveWcov = false;     % Set to false to save memory when saving final results
segLength = Inf;      % Optional speedup to cut trials into smaller segments during fitting

%% ===================
% 1c) Fit smDLAG model
% ====================

% Optional speedup to cut trials into smaller segments during fitting
seqTrueCut = cutTrials(seqTrue, 'segLength', segLength);

[estParams,~,trackedParams,flags] ...
    = em_smdlag(initParams, ...
                seqTrueCut, ...
                xDim_fit, ...
                'prior', prior, ...
                'tol', tol, ...
                'maxIters', maxIters, ...
                'freqLB', freqLB, ...
                'freqParam', freqParam, ...
                'learnTau', learnTau, ...
                'learnDelays', learnDelays, ...
                'learnInducingLocs', learnInducingLocs, ...
                'verbose', verbose, ...
                'minVarFrac', minVarFrac, ...
                'maxDelayFrac', maxDelayFrac, ...
                'maxTauFrac', maxTauFrac, ...
                'pruneX', pruneX, ...
                'saveCcov', saveCcov, ...
                'saveWcov', saveWcov);
           
%% =======================
% 1d) Save fitting results
% ========================
           
save('demo/results/demo_smdlag_results.mat', ...
     'estParams', 'trackedParams', 'flags');

%% =======================
% 2a) Load fitting results
% ========================

load('demo/results/demo_smdlag_results.mat');

%% ========================
% 2b) Check fitting results
% =========================

% Display flags indicating fitting procedure status
% flags.Convergence -- Indicates that fitting converged according to 'tol'.
%                      Else 'maxIters' was reached before convergence.
% flags.DecreasingLowerBound -- Indicates that the lower bound (objective
%                               function) decreased at some point during
%                               fitting, which suggests there's something
%                               wrong with the fitting or data.
% flags.PrivateVarianceFloor -- Indicates that the variance floor was used on
%                               one or more observed dimensions. 
%                               See em_mdlag.m header for more info.
% flags.xDimsRemoved         -- Number of latent dimensions removed (if 
%                               pruneX is true) due to low variance in all 
%                               groups.
fprintf('\n');
disp(flags)
plotFittingProgress(trackedParams, binWidth, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

%% ==============================================================
% 3) Visualize recovery of the loading matrix and ARD parameters
% ==============================================================

% Loadings matrices

% Ground truth
Ctrue = vertcat(paramsTrue.Cs{:});
hinton(Ctrue);

% Estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue.
Cest = vertcat(estParams.C.means{:});
hinton(Cest);

% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.

% Ground truth
alpha_inv_true = 1./paramsTrue.alphas;
% Normalize by the shared variance in each area
alpha_inv_rel_true = alpha_inv_true ./ sum(alpha_inv_true,2);
hinton(alpha_inv_rel_true);

% Estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not match
%       the order in alpha_true.
alpha_inv_est = 1./estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);
hinton(alpha_inv_rel_est);

%% ========================================================================
% 4) Explore leave-group-out prediction performance and SNRs on train data
% ========================================================================

% Leave-group-out predictive performance
[R2, MSE] = pred_smdlag(seqTrue, estParams)

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
snr = computeSNR(estParams.C, estParams.phi)

%% =================================================================
% 5) Determine and then visualize the dimensionalities of all types
% =================================================================

cutoff_sharedvar = 0.01; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be signficiant
[dims,sigDims,varExp,dimTypes] = computeDimensionalities(estParams, ...
                                                         cutoff_sharedvar, ...
                                                         cutoff_snr);
                                              
% Visualize the number of each type of dimension
plotDimensionalities(dims, dimTypes, ...
                     'groupNames', groupNames, ...
                     'plotZeroDim', false);

% Visualize the shared variance explained by each dimension type in each
% group
plotVarExp(varExp, dimTypes, ...
           'groupNames', groupNames, ...
           'plotZeroDim', false);

%% ======================================
% 6) Visualize recovery of GP parameters
% ======================================

% Ground truth
plotGPparams_smdlag(paramsTrue,binWidth,'sigDims',sigDimsTrue,'units',units);

% Estimate
plotGPparams_smdlag(estParams,binWidth,'sigDims',sigDims,'units',units);

%% ============================================
% 7) Visualize recovery of latent time courses
% ============================================

xspec = 'xsm'; % 'xve' gives latent time courses scales by shared variance
               % 'xsm' gives latent time courses with normalized variances
                 
% Ground truth
[seqTrue, sortParams] = scaleByVarExp(seqTrue, paramsTrue, alpha_inv_rel_true, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_smdlag(seqTrue, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {}, ...
                     'plotW', true);
                 
% Estimate
[seqEst,~,~] = inferXW(seqTrue, estParams);
% NOTE: In general, latent time courses are unordered, and will not match
%       the order in seqTrue.
[seqEst, sortParams] = scaleByVarExp(seqEst, estParams, alpha_inv_rel_est, ...
                                     'sortDims', false, ...
                                     'sortGroup', 1, ...
                                     'numDim', 10, ...
                                     'indatafield', 'xsm', ...
                                     'outdatafield', 'xve');
plotDimsVsTime_smdlag(seqEst, xspec, sortParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {}, ...
                     'trialColors', {}, ...
                     'plotW', true);

%% =============================================================
% 8) Perform a pairwise analysis of interactions between groups
% =============================================================

% Compute pairwise dimensionalities and shared variances explained
[pairDims,pairVarExp,pairs] = computeDims_pairs(dims,dimTypes,varExp);

% Visualize results
plotDims_pairs(pairDims, pairs, numGroups, 'groupNames', groupNames);
plotVarExp_pairs(pairVarExp, pairs, numGroups, 'groupNames', groupNames);
