% results_summary_basic_demo.m
%
% Description: This script compares mDLAG-time, mDLAG-inducing, and
%              mDLAG-frequency results on the basic demo data.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_basic_demo;

%% Load fitted model and ground truth data

% Ground truth data
currDataFile = sprintf('%s/%s', dataDir, dataFile);
load(dataFile);

% Time domain fit
resultFile_time = sprintf('%s/mdlag_results_basic_demo_time.mat', resultDir);
results_time = load(resultFile_time);

% Inducing point fit
resultFile_sparse = sprintf('%s/smdlag_S050_results_basic_demo.mat', resultDir);
results_sparse = load(resultFile_sparse);
% Restructure the smDLAG structures to be compatible with the others
results_sparse.trackedParams.gp_params.Ds = results_sparse.trackedParams.Ds;
results_sparse.trackedParams.gp_params.gams = results_sparse.trackedParams.gams;
results_sparse.trackedParams = rmfield(results_sparse.trackedParams, {'Ds', 'gams'});

% Frequency domain fit
resultFile_freq = sprintf('%s/mdlag_results_basic_demo_freq.mat', resultDir);
results_freq = load(resultFile_freq);

%% Inspect the fitting progress of each method

% Time domain
plotFittingProgress(results_time.trackedParams, ...
                    binWidth, ...
                    results_time.estParams.covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);

% Inducing points
plotFittingProgress(results_sparse.trackedParams, ...
                    binWidth, ...
                    results_sparse.estParams.covType, ...
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

%% Compare lower bound progress over clock time

figure;
hold on;
% Time domain
h1 = plot(cumsum(results_time.trackedParams.iterTime), ...
          results_time.trackedParams.lb, ...
          'color', TIMECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% Inducing points
h2 = plot(cumsum(results_sparse.trackedParams.iterTime), ...
          results_sparse.trackedParams.lb, ...
          'color', SPARSECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% Frequency domain
h3 = plot(cumsum(results_freq.trackedParams.iterTime), ...
          results_freq.trackedParams.lb, ...
          'color', FREQCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
xscale('log');
xlabel('Elapsed clock time (s)');
ylabel('Lower bound value');
axis square;
legend([h1 h2 h3], 'Time domain', 'Inducing points', 'Freq. domain', ...
    'location', 'southeast');

% Runtime summaries
fprintf('Avg. runtime per iteration (s):\n')
fprintf('    time domain:    %f\n', mean(results_time.trackedParams.iterTime))
fprintf('    inducing:       %f\n', mean(results_sparse.trackedParams.iterTime))
fprintf('    freq domain:    %f\n', mean(results_freq.trackedParams.iterTime))

fprintf('Total iterations:\n')
fprintf('    time domain:    %d\n', length(results_time.trackedParams.iterTime))
fprintf('    inducing:       %d\n', length(results_sparse.trackedParams.iterTime))
fprintf('    freq domain:    %d\n', length(results_freq.trackedParams.iterTime))

fprintf('Total runtime (s):\n')
fprintf('    time domain:    %f\n', sum(results_time.trackedParams.iterTime))
fprintf('    inducing:       %f\n', sum(results_sparse.trackedParams.iterTime))
fprintf('    freq domain:    %f\n', sum(results_freq.trackedParams.iterTime))

%% Visualize recovery of the loading matrix and ARD parameters

% Loadings matrices

% Ground truth
Ctrue = vertcat(trueParams.Cs{:});
hinton(Ctrue);

% Time domain estimate
% NOTE: In general, the columns of Cest are unordered, and will not match
%       the order in Ctrue. Here, we can reorder the latents to 
%       facilitate comparison with the ground truth.
Cest = vertcat(results_time.estParams.C.means{:});
reorder_time = [1  2  7  4];
rescale_time = [1  1 -1  1];
hinton(Cest(:,reorder_time).*rescale_time);
% Compute normalized distance from the ground truth
Cerr_time = norm(Ctrue - Cest(:,reorder_time).*rescale_time, 'fro') ...
    / norm(Ctrue, 'fro');

% Inducing point estimate
Cest = vertcat(results_sparse.estParams.C.means{:});
reorder_sparse = [1  2  7  5];
rescale_sparse = [1  1 -1  1];
hinton(Cest(:,reorder_sparse).*rescale_sparse);
% Compute distance from the ground truth
Cerr_sparse = norm(Ctrue - Cest(:,reorder_sparse).*rescale_sparse, 'fro') ...
    / norm(Ctrue, 'fro');

% Frequency domain estimate
Cest = vertcat(results_freq.estParams.C.means{:});
reorder_freq = [1  2  7  5];
rescale_freq = [1  1 -1  1];
hinton(Cest(:,reorder_freq).*rescale_freq);
% Compute distance from the ground truth
Cerr_freq = norm(Ctrue - Cest(:,reorder_freq).*rescale_freq, 'fro') ...
    / norm(Ctrue, 'fro');

% Report performance summaries
fprintf('Normalized loading matrix estimation error:\n')
fprintf('    time domain:    %0.5f\n', Cerr_time)
fprintf('    inducing:       %0.5f\n', Cerr_sparse)
fprintf('    freq domain:    %0.5f\n', Cerr_freq)

% Alpha parameters
% The following plots visualize the shared variance explained by each
% latent variable in each area.

% Ground truth
alpha_inv_true = 1./trueParams.alphas;
% Normalize by the shared variance in each area
alpha_inv_rel_true = alpha_inv_true ./ sum(alpha_inv_true,2);

figure;
hold on;
b = bar(alpha_inv_rel_true');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.6]);
title('Ground truth')

% Time domain estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not match
%       the order in alpha_true.
alpha_inv_est = 1./results_time.estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

figure;
hold on;
b = bar(alpha_inv_rel_est(:,reorder_time)');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.6]);
title('Time domain estimate')

% Inducing point estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not match
%       the order in alpha_true.
alpha_inv_est = 1./results_sparse.estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

figure;
hold on;
b = bar(alpha_inv_rel_est(:,reorder_sparse)');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.6]);
title('Inducing point estimate')

% Frequency domain estimate
% NOTE: In general, the columns of alpha_est are unordered, and will not match
%       the order in alpha_true.
alpha_inv_est = 1./results_freq.estParams.alpha.mean;
% Normalize by the shared variance in each area
alpha_inv_rel_est = alpha_inv_est ./ sum(alpha_inv_est,2);

figure;
hold on;
b = bar(alpha_inv_rel_est(:,reorder_freq)');
for groupIdx = 1:numGroups
    b(groupIdx).FaceColor = groupColors{groupIdx};
end
xlabel('Latent variable');
ylabel('Frac. shared var. exp.');
ylim([0 0.6]);
title('Freq. domain estimate')

%% Visualize recovery of GP parameters

% Ground truth
gp_params_true = getGPparams_mdlag(trueParams, binWidth);

% Time domain estimate
gp_params_time = getGPparams_mdlag(results_time.estParams, binWidth);

% Inducing points estimate
gp_params_sparse = getGPparams_mdlag(results_sparse.estParams, binWidth);

% Frequency domain estimate
gp_params_freq = getGPparams_mdlag(results_freq.estParams, binWidth);

% Display timescales
figure;
hold on;
b = bar([gp_params_true.tau; ...
         gp_params_time.tau(reorder_time); ...
         gp_params_sparse.tau(reorder_freq); ...
         gp_params_freq.tau(reorder_freq)]');
b(1).FaceColor = GTCOLOR;
b(2).FaceColor = TIMECOLOR;
b(3).FaceColor = SPARSECOLOR;
b(4).FaceColor = FREQCOLOR;
xlabel('Latent variable');
ylabel('Timescale (ms)');
ylim([0 max(gp_params_true.tau)+5]);
xlim([0.5, xDim+0.5]);
xticks(1:xDim);
axis square;
legend([b(1), b(2), b(3), b(4)], {'Ground truth', 'Time est.', ...
    'Inducing est.', 'Freq. est.'}, 'location', 'southoutside');

% Display time delays
figure;
hold on;
b = bar([gp_params_true.D(2,1:2); ...
         gp_params_time.D(2,reorder_time(1:2)); ...
         gp_params_sparse.D(2,reorder_sparse(1:2)); ...
         gp_params_freq.D(2,reorder_freq(1:2))]');
b(1).FaceColor = GTCOLOR;
b(2).FaceColor = TIMECOLOR;
b(3).FaceColor = SPARSECOLOR;
b(4).FaceColor = FREQCOLOR;
xlabel('Latent variable');
ylabel('Time delay (ms)');
ylim([-25 25]);
xlim([0.5, 2.5]);
xticks(1:2);
axis square;
legend([b(1), b(2), b(3), b(4)], {'Ground truth', 'Time est.', ...
    'Inducing est.', 'Freq. est.'}, 'location', 'southoutside');

%% Visualize recovery of latent time courses

trialIdx = 1;  % Choose an example trial
                        
% Time domain estimates
[seqEst_time,~,~] = inferX(seqTrue, results_time.estParams);

% Inducing point estimates
[seqEst_sparse,~,~] = inferX(seqTrue, results_sparse.estParams);

% Frequency domain estimates
[seqEst_freq,~,~] = inferX(seqTrue, results_freq.estParams);

% Overlay estimates on ground truth, for an example trial                 
seqEst_time(trialIdx).xsm ...
    = seqEst_time(trialIdx).xsm([reorder_time reorder_time+results_time.estParams.xDim],:) .* repmat(rescale_time, 1, 2)'; 
seqEst_sparse(trialIdx).xsm ...
    = seqEst_sparse(trialIdx).xsm([reorder_sparse reorder_sparse+results_sparse.estParams.xDim],:) .* repmat(rescale_sparse, 1, 2)';
seqEst_freq(trialIdx).xsm ...
    = seqEst_freq(trialIdx).xsm([reorder_freq reorder_freq+results_freq.estParams.xDim],:) .* repmat(rescale_freq, 1, 2)';
plotDimsVsTime_mdlag([seqTrue(trialIdx); seqEst_time(trialIdx); seqEst_freq(trialIdx); seqEst_sparse(trialIdx)], ...
                     'xsm', ...
                     trueParams, ...
                     binWidth, ...
                     'nPlotMax', 1, ...
                     'plotSingle', false, ...
                     'plotMean', true, ...
                     'units', units, ...
                     'trialGroups', {[1], [2], [3], [4]}, ...
                     'trialColors', {GTCOLOR, TIMECOLOR, FREQCOLOR, SPARSECOLOR});

%% Compare leave-group-out prediction performance and SNRs on train data

fprintf('Leave-group-out R^2:\n');

% Time domain performance
[R2_time, ~] = pred_mdlag(seqTrue, results_time.estParams);
fprintf('    time domain:    %1.4f\n', R2_time);

% Inducing point performance
[R2_sparse, ~] = pred_mdlag(seqTrue, results_sparse.estParams);
fprintf('    inducing:       %1.4f\n', R2_sparse);

% Frequency domain performance
[R2_freq, ~] = pred_mdlag(seqTrue, results_freq.estParams);
fprintf('    freq domain:    %1.4f\n', R2_freq);

% Signal-to-noise ratio of each group, according to estimated mDLAG model
% parameters
fprintf('Signal-to-noise ratio of each group:\n');

% Time domain
snr = computeSNR(results_time.estParams.C, results_time.estParams.phi);
fprintf('    time domain:    %1.4f  %1.4f  \n', snr);

% Inducing points
snr = computeSNR(results_sparse.estParams.C, results_sparse.estParams.phi);
fprintf('    inducing:       %1.4f  %1.4f  \n', snr);

% Frequency domain
snr = computeSNR(results_freq.estParams.C, results_freq.estParams.phi);
fprintf('    freq domain:    %1.4f  %1.4f  \n', snr);
                 
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

% Inducing points
[seqEst_sparse,~,~] = inferX(seqTrue, results_sparse.estParams);

Xs_gt = seq2cell2D(seqTrue, trueParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');
Xs_est = seq2cell2D(seqEst_sparse, results_sparse.estParams.xDim.*ones(1,numGroups), 'datafield', 'xsm');

Y_gt = blkdiag(trueParams.Cs{:}) * vertcat(Xs_gt{:}) + vertcat(trueParams.ds{:});
Y_est = blkdiag(results_sparse.estParams.C.means{:}) * vertcat(Xs_est{:}) + results_sparse.estParams.d.mean;
RSS = sum( sum( ( Y_gt - Y_est ).^2, 1 ) );
TSS = sum( sum( ( Y_gt - repmat( mean(Y_gt,2), [1 size(Y_gt,2)] ) ).^2, 1 ) );
R2_sparse = 1 - RSS / TSS;
fprintf('    inducing:       %1.4f\n', R2_sparse);

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
