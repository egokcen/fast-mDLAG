% ===================================================================
% mDLAG-frequency DEMO: Run this script from the main mDLAG directory
% ===================================================================
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% ================
% 0) Load demo data
% =================

% Synthetic data generated from a mDLAG model.
% See simulation/mdlag_groundtruth_data.m for a demonstration of how to
% generate simulated data from a mDLAG model.
dat_file = 'demo/data/mdlag-freq_demo_data';
fprintf('Reading from %s \n',dat_file);
load(dat_file);
numGroups = length(paramsTrue.yDims);
yDims = paramsTrue.yDims;
units = 'ms';

% Color scheme
GTCOLOR = '#7F7F7F';
TIMECOLOR = 'k';
FREQCOLOR = '#D35FBC';

%% ====================================
% 1a) Initialize mDLAG model parameters
% =====================================

xDim_fit = 2;            % Number of latent variables to be fitted
covType = 'rbf';         % Type of GP covariance kernel
randomSeed = 0;          % Seed the random number generator for reproducibility
prior_val = 1e-12;       % Set to very small value to set uninformative prior
prior.d.beta = prior_val;
prior.phi.a = prior_val;
prior.phi.b = prior_val;
prior.alpha.a = prior_val;
prior.alpha.b = prior_val;
saveCcov = false;        % Set to false to save memory when saving initParams
initParams = init_mdlag(seqTrue, ...
                        yDims, ...
                        xDim_fit, ...
                        binWidth, ...
                        'covType', covType, ...
                        'prior', prior, ...
                        'randomSeed', randomSeed, ...
                        'saveCcov', saveCcov);

%% ===================================================================
% 1b) Compare the runtimes of time domain and frequency domain fitting
% ====================================================================

% Common arguments across both approaches

% Let's explicitly define relevant optional arguments, for 
% the sake of demonstration:
tol = 1e-8;           % Tolerance to determine EM convergence
maxIters = 10;        % Maximum EM iterations
freqLB = 1;           % Check for convergence of lower bound every freqLB EM iterations
freqParam = 10;       % Store intermediate delay and timescale estimates every freqParam EM iterations
learnDelays = true;   % Toggle whether to learn delay parameters
verbose = true;       % Print fitting progress
minVarFrac = 0.001;   % Private noise variances will not be allowed to go below this value
maxDelayFrac = 0.5;   % Maximum delay magnitude (unit: fraction of trial length)
maxTauFrac = 1.0;     % Maximum timescale magnitude (unit: fraction of trial length)
pruneX = true;        % For speed-up, remove latents that become inactive in all groups
saveXcov = false;     % Set to false to save memory when saving final results
saveCcov = false;     % Set to false to save memory when saving final results

% As a preprocessing step, compute the FFT of the observed data
seqTrue = fftseq(seqTrue, 'y', 'yfft');

%% =====================================================================
% 1c) Time domain fitting
%     NOTE: You can skip this section if a trained model already exists.
% ======================================================================

[estParams,~,trackedParams,flags] ...
    = em_mdlag(initParams, ...
               seqTrue, ...
               xDim_fit, ...
               'prior', prior, ...
               'tol', tol, ...
               'maxIters', maxIters, ...
               'freqLB', freqLB, ...
               'freqParam', freqParam, ...
               'learnDelays', learnDelays, ...
               'verbose', verbose, ...
               'minVarFrac', minVarFrac, ...
               'maxDelayFrac', maxDelayFrac, ...
               'maxTauFrac', maxTauFrac, ...
               'pruneX', pruneX, ...
               'saveXcov', saveXcov, ...
               'saveCcov', saveCcov);

% Save fitting results
save('demo/results/demo_mdlag_freq_runtime_time.mat', ...
     'estParams', 'trackedParams', 'flags');

%% =====================================================================
% 1d) Frequency domain fitting
%     NOTE: You can skip this section if a trained model already exists.
% ======================================================================

[estParams,~,trackedParams,flags] ...
    = em_mdlag_freq(initParams, ...
                    seqTrue, ...
                    xDim_fit, ...
                    'prior', prior, ...
                    'tol', tol, ...
                    'maxIters', maxIters, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'learnDelays', learnDelays, ...
                    'verbose', verbose, ...
                    'minVarFrac', minVarFrac, ...
                    'maxDelayFrac', maxDelayFrac, ...
                    'maxTauFrac', maxTauFrac, ...
                    'pruneX', pruneX, ...
                    'saveXcov', saveXcov, ...
                    'saveCcov', saveCcov);

% Save fitting results
save('demo/results/demo_mdlag_freq_runtime_freq.mat', ...
     'estParams', 'trackedParams', 'flags');

%% =======================
% 1e) Investigate runtimes
% ========================

fprintf('\nAvg. runtime per iteration (s):\n')

% Time domain
load('demo/results/demo_mdlag_freq_runtime_time.mat');
fprintf('    time domain:    %f\n', mean(trackedParams.iterTime))

% Frequency domain
load('demo/results/demo_mdlag_freq_runtime_freq.mat');
fprintf('    freq domain:    %f\n', mean(trackedParams.iterTime))

%% ======================================================================
% 2a) Fully fit a mDLAG model using the frequency domain approach
%     NOTE: You can skip to Section 2c if a trained model already exists.
% =======================================================================

maxIters = 5000;

[estParams,~,trackedParams,flags] ...
    = em_mdlag_freq(initParams, ...
                    seqTrue, ...
                    xDim_fit, ...
                    'prior', prior, ...
                    'tol', tol, ...
                    'maxIters', maxIters, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'learnDelays', learnDelays, ...
                    'verbose', verbose, ...
                    'minVarFrac', minVarFrac, ...
                    'maxDelayFrac', maxDelayFrac, ...
                    'maxTauFrac', maxTauFrac, ...
                    'pruneX', pruneX, ...
                    'saveXcov', saveXcov, ...
                    'saveCcov', saveCcov);

%% =======================
% 2b) Save fitting results
% ========================

save('demo/results/demo_mdlag_freq_fullfit.mat', ...
     'estParams', 'trackedParams', 'flags');

%% =======================
% 2c) Load fitting results
% ========================

load('demo/results/demo_mdlag_freq_fullfit.mat');

%% ========================
% 2d) Check fitting results
% =========================
% 
% Display flags indicating fitting procedure status
% flags.Convergence -- Indicates that fitting converged according to 'tol'.
%                      Else 'maxIters' was reached before convergence.
% flags.DecreasingLowerBound -- Indicates that the lower bound (objective
%                               function) decreased at some point during
%                               fitting, which suggests there's something
%                               wrong with the fitting or data.
% flags.PrivateVarianceFloor -- Indicates that the variance floor was used
%                               on one or more observed dimensions. 
%                               See em_mdlag.m header for more info.
% flags.xDimsRemoved         -- Number of latent dimensions removed (if 
%                               pruneX is true) due to low variance in all 
%                               groups.
fprintf('\n');
disp(flags)
plotFittingProgress(trackedParams, binWidth, estParams.covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

%% ======================================
% 2e) Visualize recovery of GP parameters
% =======================================

% Ground truth
plotGPparams_mdlag(paramsTrue,binWidth,'units',units);

% Estimate
plotGPparams_mdlag(estParams,binWidth,'units',units);

%% ============================================
% 2f) Visualize recovery of latent time courses
% =============================================

% NOTE: We fit a mDLAG model via the frequency domain. Here, we'll use
%       those fitted parameters, but perform all inference in the time
%       domain.

trialIdx = 1;  % Choose an example trial to plot

% Ground truth
plotDimsVsTime_mdlag(seqTrue(trialIdx), 'xsm', paramsTrue, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1]}, ...
                     'trialColors', {GTCOLOR});
                 
% Estimate
[seqEst,~,~] = inferX(seqTrue, estParams);
% NOTE: In general, latent time courses are unordered, and will not match
%       the order in seqTrue.
plotDimsVsTime_mdlag(seqEst(trialIdx), 'xsm', estParams, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1]}, ...
                     'trialColors', {FREQCOLOR});
                 
% Overlay estimate on ground truth, for an example trial  
reorder = [2 1];    % Reorder estimated latents as needed
rescale = [-1 -1];  % Flip estimated latents as needed
seqEst(trialIdx).xsm = seqEst(trialIdx).xsm([reorder reorder+estParams.xDim],:) .* repmat(rescale, 1, 2)';                
plotDimsVsTime_mdlag([seqTrue(trialIdx); seqEst(trialIdx)], 'xsm', paramsTrue, binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1], [2]}, ...
                     'trialColors', {GTCOLOR, FREQCOLOR});

%% ======================================================================
% 3a) Using the same set of parameters, compare time and frequency domain
%     inference
% =======================================================================

trialIdx = 1;  % Choose an example trial to plot

% Load fitted parameters from Section 2
load('demo/results/demo_mdlag_freq_fullfit.mat');

fprintf('\nLatent inference runtime (s):\n');

% Time domain estimate
tic;
[seqEst_time,~,~] = inferX(seqTrue, estParams);
runtime = toc;
fprintf('    time domain:  %1.4f\n', runtime);

% Frequency domain estimate
% Take the FFT of each neuron's activity on each trial
seqTrue = fftseq(seqTrue,'y','yfft');
% Perform inference in the frequency domain
tic;
[seqEst_freq,~,~] = inferX_freq(seqTrue, estParams);
runtime = toc;
fprintf('    freq_domain:  %1.4f\n', runtime);
% Convert estimated latents back to the time domain
seqEst_freq = freq2time_mdlag(seqEst_freq, estParams, ...
    'infield', 'xfft', 'outfield', 'xsm');

% Let's plot ground truth, time and frequency domain estimates on top of 
% each other
latentIdx = 2;      % Choose one latent to plot
reorder = [2 1];    % Reorder estimated latents as needed
rescale = [-1 -1];  % Flip estimated latents as needed

% Notice the edge effects that appear in the frequency domain estimates
figure;
hold on;
% Ground truth
h1 = plot(1:length(seqTrue(trialIdx).xsm(latentIdx,:)), ...
          seqTrue(trialIdx).xsm(latentIdx,:), ...
          'color', GTCOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
% Time domain approach
h2 = plot(1:length(seqEst_time(trialIdx).xsm(reorder(latentIdx),:)), ...
          seqEst_time(trialIdx).xsm(reorder(latentIdx),:)*rescale(latentIdx), ...
          'color', TIMECOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
% Frequency domain approach
h3 = plot(1:length(seqEst_freq(trialIdx).xsm(reorder(latentIdx),:)), ...
          seqEst_freq(trialIdx).xsm(reorder(latentIdx),:)*rescale(latentIdx), ...
          'color', FREQCOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
legend([h1, h2, h3], {'Ground truth', 'Time', 'Freq'});
xlabel('Time point');
ylabel('x');

% Error between time and frequency domain approaches as a function of time.
% Again, note the error increase at the edges.
Xtime = cat(3,seqEst_time.xsm);
Xfreq = cat(3,seqEst_freq.xsm);
Xerr = sqrt(mean((Xtime(reorder,:,:) - Xfreq(reorder,:,:)).^2,3)); % RMSE
figure;
for j = 1:size(Xerr,1)
    subplot(size(Xerr,1),1,j);
    hold on;
    plot(Xerr(j,:),'k-', 'linewidth', 1.5);
    xlabel('Time');
    ylabel('Error');
    title(sprintf('Latent %d', j));
end

%% ======================================================================
% 3b) Using the same set of parameters, compare time and frequency domain
%     prediction
% =======================================================================

% Load fitted parameters from Section 2
load('demo/results/demo_mdlag_freq_fullfit.mat');

% Leave-group-out predictive performance

fprintf('\nLeave-group-out prediction:\n')

% Time domain
tic;
[R2, ~] = pred_mdlag(seqTrue, estParams);
runtime = toc;
fprintf('Using time domain:\n');
fprintf('    R^2:          %1.4f\n', R2);
fprintf('    runtime (s):  %1.4f\n', runtime);

% Frequency domain
tic;
[R2, ~] = pred_mdlag_freq(seqTrue, estParams);
runtime = toc;
fprintf('Using freq domain:\n');
fprintf('    R^2:          %1.4f\n', R2);
fprintf('    runtime (s):  %1.4f\n', runtime);