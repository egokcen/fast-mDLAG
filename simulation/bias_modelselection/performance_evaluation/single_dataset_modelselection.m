% single_dataset_modelselection.m
%
% Description: Inspect a single dataset in detail.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_modelselection;

%% Load fitted model and ground truth data

runIdx = 30;
T = Tlist(4);

% Ground truth data
currDataFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
load(currDataFile);
% Time domain fit
results_time = load( ...
    sprintf('%s/run%03d/T%04d_time.mat', resultDir, runIdx, T) ...
);
% Frequency domain fit
results_freq = load( ...
    sprintf('%s/run%03d/T%04d_freq.mat', resultDir, runIdx, T) ...
);
% Frequency domain fit with taper
results_taper = load( ...
    sprintf('%s/run%03d/T%04d_freq_hamming.mat', resultDir, runIdx, T) ...
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

%% Compare SNRs on train data

% Signal-to-noise ratio, according to estimated mDLAG model parameters
fprintf('Signal-to-noise ratio of each group:\n');

% Time domain
snr = computeSNR(results_time.estParams.C, results_time.estParams.phi);
fprintf('    time domain:    %1.4f \n', snr);

% Frequency domain
snr = computeSNR(results_freq.estParams.C, results_freq.estParams.phi);
fprintf('    freq domain:    %1.4f  \n', snr);

% Frequency domain with taper
snr = computeSNR(results_taper.estParams.C, results_taper.estParams.phi);
fprintf('    freq + taper:   %1.4f \n', snr);
        
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
