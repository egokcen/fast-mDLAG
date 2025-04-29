% single_dataset_samplingrate.m
%
% Description: Inspect a single dataset in detail.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_bias_samplingrate;

%% Load fitted model and ground truth data

runIdx = 1;
binWidth = binWidthList(end);

% Ground truth data
currDataFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
load(currDataFile);
% Time domain fit
results_time = load( ...
    sprintf('%s/run%03d/period%03d_time.mat', resultDir, runIdx, binWidth) ...
);
% Frequency domain fit
results_freq = load( ...
    sprintf('%s/run%03d/period%03d_freq.mat', resultDir, runIdx, binWidth) ...
);
% Frequency domain fit with taper
results_taper = load( ...
    sprintf('%s/run%03d/period%03d_freq_hamming.mat', resultDir, runIdx, binWidth) ...
);

%% Inspect the fitting progress of each method

% Time domain
plotFittingProgress(results_time.trackedParams, ...
                    binWidth, ...
                    results_time.estParams.covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

% Frequency domain
plotFittingProgress(results_freq.trackedParams, ...
                    binWidth, ...
                    results_freq.estParams.covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

% Frequency domain with taper
plotFittingProgress(results_taper.trackedParams, ...
                    binWidth, ...
                    results_freq.estParams.covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

%% Compare lower bound progress over clock time

figure;
hold on;
% Time domain
h1 = plot(cumsum(results_time.trackedParams.iterTime), ...
          results_time.trackedParams.lb, ...
          'color', TIMECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% Frequency domain
h2 = plot(cumsum(results_freq.trackedParams.iterTime), ...
          results_freq.trackedParams.lb, ...
          'color', FREQCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% Frequency domain with taper
h3 = plot(cumsum(results_taper.trackedParams.iterTime), ...
          results_taper.trackedParams.lb, ...
          'color', TAPERCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
xscale('log');
xlabel('Elapsed clock time (s)');
ylabel('Lower bound value');
legend([h1 h2 h3], 'Time domain', 'Freq. domain', 'Freq. domain + taper', ...
    'location', 'southeast');

% Runtime summaries
fprintf('Avg. runtime per iteration (s):\n')
fprintf('    time domain:    %f\n', mean(results_time.trackedParams.iterTime))
fprintf('    freq domain:    %f\n', mean(results_freq.trackedParams.iterTime))
fprintf('    freq + taper:   %f\n', mean(results_taper.trackedParams.iterTime))

fprintf('Total iterations:\n')
fprintf('    time domain:    %d\n', length(results_time.trackedParams.iterTime))
fprintf('    freq domain:    %d\n', length(results_freq.trackedParams.iterTime))
fprintf('    freq + taper:   %d\n', length(results_taper.trackedParams.iterTime))

fprintf('Total runtime (s):\n')
fprintf('    time domain:    %f\n', sum(results_time.trackedParams.iterTime))
fprintf('    freq domain:    %f\n', sum(results_freq.trackedParams.iterTime))
fprintf('    freq + taper:   %f\n', sum(results_taper.trackedParams.iterTime))

%% Visualize recovery of the loading matrix

% Loadings matrices

% Ground truth
Ctrue = vertcat(trueParams.Cs{:});
hinton(Ctrue);

% Time domain estimate
rescale_time = -1;
Cest = vertcat(results_time.estParams.C.means{:});
hinton(Cest.*rescale_time);

% Frequency domain estimate
rescale_freq = -1;
Cest = vertcat(results_freq.estParams.C.means{:});
hinton(Cest.*rescale_freq);

% Frequency domain + taper
rescale_taper = -1;
Cest = vertcat(results_taper.estParams.C.means{:});
hinton(Cest.*rescale_taper);

%% Visualize recovery of latent time courses

trialIdx = 1;  % Choose an example trial
                        
% Time domain estimates
[seqEst_time,~,~] = inferX(seqTrue, results_time.estParams);
seqEst_time(trialIdx).xsm = seqEst_time(trialIdx).xsm .* rescale_time;

% Frequency domain estimates
[seqEst_freq,~,~] = inferX(seqTrue, results_freq.estParams);
seqEst_freq(trialIdx).xsm = seqEst_freq(trialIdx).xsm .* rescale_freq;

% Frequency domain with taper
[seqEst_taper,~,~] = inferX(seqTrue, results_taper.estParams);
seqEst_taper(trialIdx).xsm = seqEst_taper(trialIdx).xsm .* rescale_taper;

% Overlay estimates on ground truth, for an example trial                 
plotDimsVsTime_mdlag([seqTrue(trialIdx); seqEst_time(trialIdx); seqEst_freq(trialIdx); seqEst_taper(trialIdx)], ...
                     'xsm', ...
                     trueParams, ...
                     min(binWidthList), ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1], [2], [3], [4]}, ...
                     'trialColors', {GTCOLOR, TIMECOLOR, FREQCOLOR, TAPERCOLOR});

%% Compare leave-group-out prediction performance and SNRs on train data

fprintf('Leave-group-out R^2:\n');

% Time domain performance
[R2_time, ~] = pred_mdlag(seqTrue, results_time.estParams);
fprintf('    time domain:    %1.4f\n', R2_time);

% Frequency domain performance
[R2_freq, ~] = pred_mdlag(seqTrue, results_freq.estParams);
fprintf('    freq domain:    %1.4f\n', R2_freq);

% Frequency domain with taper
[R2_taper, ~] = pred_mdlag(seqTrue, results_taper.estParams);
fprintf('    freq + taper:   %1.4f\n', R2_taper);

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
fprintf('Signal-to-noise ratio of each group:\n');

% Time domain
snr = computeSNR(results_time.estParams.C, results_time.estParams.phi);
fprintf('    time domain:    %1.4f  %1.4f  \n', snr);

% Frequency domain
snr = computeSNR(results_freq.estParams.C, results_freq.estParams.phi);
fprintf('    freq domain:    %1.4f  %1.4f  \n', snr);

% Frequency domain with taper
snr = computeSNR(results_taper.estParams.C, results_taper.estParams.phi);
fprintf('    freq + taper:   %1.4f  %1.4f  \n', snr);
        
%% Compare recovery of latent time courses

fprintf('R^2, latent reconstruction:\n');

% Time domain
[seqEst_time,~,~] = inferX(seqTrue, results_time.estParams);

Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst_time, results_time.estParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = blkdiag(results_time.estParams.C.means{:}) * vertcat(Xs_est{:}) + results_time.estParams.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_time = 1 - RSS / TSS;
fprintf('    time domain:    %1.4f\n', R2_time);

% Frequency domain
[seqEst_freq,~,~] = inferX(seqTrue, results_freq.estParams);

Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst_freq, results_freq.estParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = blkdiag(results_freq.estParams.C.means{:}) * vertcat(Xs_est{:}) + results_freq.estParams.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_freq = 1 - RSS / TSS;
fprintf('    freq domain:    %1.4f\n', R2_freq);

% Frequency domain with taper
[seqEst_taper,~,~] = inferX(seqTrue, results_taper.estParams);

Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst_taper, results_taper.estParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = blkdiag(results_taper.estParams.C.means{:}) * vertcat(Xs_est{:}) + results_taper.estParams.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_taper = 1 - RSS / TSS;
fprintf('    freq + taper:   %1.4f\n', R2_taper);
