% bias_modelselection.m
%
% Description: Characterize the bias of mDLAG dimensionality estimates as a 
%              function of sequence (trial) length. Compare time, inducing
%              point, and frequency domain fitting approaches.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_modelselection;

%% Characterize estimates of dimensionality
cutoff_sharedvar = 0.01; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be signficiant
SIG_IDX = 2;             % Access significant dimensions for one group

% Time domain
dim_time = nan(numPartitions, snrPartitions, runsPerSNR);

% Inducing points
dim_sparse = nan(numPartitions, snrPartitions, runsPerSNR);

% Frequency domain
dim_freq = nan(numPartitions, snrPartitions, runsPerSNR);

for snrIdx = 1:snrPartitions
 
    for runIdx = 1:runsPerSNR

        fprintf('Run %d of %d\n', runIdx+(snrIdx-1)*runsPerSNR, numRuns);
    
        for i = 1:numPartitions
            T = Tlist(i);

            % Load fitted models
            results_time = load( ...
                sprintf('%s/run%03d/T%04d_time.mat', ...
                        resultDir, ...
                        runIdx+(snrIdx-1)*runsPerSNR, ...
                        T) ...
            );
            results_sparse = load( ...
                sprintf('%s/run%03d/T%04d_sparse.mat', ...
                        resultDir, ...
                        runIdx+(snrIdx-1)*runsPerSNR, ...
                        T) ...
            );
            results_freq = load( ...
                sprintf('%s/run%03d/T%04d_freq.mat', ...
                        resultDir, ...
                        runIdx+(snrIdx-1)*runsPerSNR, ...
                        T) ...
            );

            % Estimated dimensionalities
            % Time domain
            [dims,~,~,~] ...
                = computeDimensionalities(results_time.estParams, ...
                                          cutoff_sharedvar, ...
                                          cutoff_snr);
            dim_time(i,snrIdx,runIdx) = dims(SIG_IDX);
            % Inducing points
            [dims,~,~,~] ...
                = computeDimensionalities(results_sparse.estParams, ...
                                          cutoff_sharedvar, ...
                                          cutoff_snr);
            dim_sparse(i,snrIdx,runIdx) = dims(SIG_IDX);
            % Frequency domain
            [dims,~,~,~] ...
                = computeDimensionalities(results_freq.estParams, ...
                                          cutoff_sharedvar, ...
                                          cutoff_snr);
            dim_freq(i,snrIdx,runIdx) = dims(SIG_IDX);

        end

    end

end

%% Visualize results

% Mean and SEM over runs
% Time domain
mean_dim_time = mean(dim_time,3);
sem_dim_time = std(dim_time,0,3) ./ sqrt(runsPerSNR);

% Inducing points
mean_dim_sparse = mean(dim_sparse,3);
sem_dim_sparse = std(dim_sparse,0,3) / sqrt(runsPerSNR);

% Frequency domain
mean_dim_freq = mean(dim_freq,3);
sem_dim_freq = std(dim_freq,0,3) / sqrt(runsPerSNR);

% Time domain
h = [];
lbls = {};
figure;
hold on;
% Ground truth
line([Tlist(1) Tlist(end)], [xDim xDim], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
% Estimate
for snrIdx = 1:snrPartitions

    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_dim_time(:,snrIdx) + sem_dim_time(:,snrIdx); ...
              flipud(mean_dim_time(:,snrIdx) - sem_dim_time(:,snrIdx))], ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = snrColors{snrIdx};
    h = [h, plot(Tlist, mean_dim_time(:,snrIdx)', ...
                 'color', snrColors{snrIdx}, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, sprintf('snr %2.2f', snrList(snrIdx))};

end
title('Time domain');
xlabel('Time points per trial');
ylabel('Dimensionality');
ylim([0, 2*xDim])
axis square;
xscale('log');
legend(h, lbls, 'location', 'southoutside');
hold off;

% Inducing points
h = [];
lbls = {};
figure;
hold on;
% Ground truth
line([Tlist(1) Tlist(end)], [xDim xDim], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
% Estimate
for snrIdx = 1:snrPartitions
 
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_dim_sparse(:,snrIdx) + sem_dim_sparse(:,snrIdx); ...
              flipud(mean_dim_sparse(:,snrIdx) - sem_dim_sparse(:,snrIdx))], ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = snrColors{snrIdx};
    h = [h, plot(Tlist, mean_dim_sparse(:,snrIdx)', ...
                 'color', snrColors{snrIdx}, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, sprintf('snr %2.2f', snrList(snrIdx))};

end
title('Inducing points');
xlabel('timep point per trial');
ylabel('Dimensionality');
ylim([0, 2*xDim])
legend(h, lbls, 'location', 'southoutside');
xscale('log');
axis square;
hold off;

% Frequency domain
h = [];
lbls = {};
figure;
hold on;
% Ground truth
line([Tlist(1) Tlist(end)], [xDim xDim], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
% Estimate
for snrIdx = 1:snrPartitions
 
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_dim_freq(:,snrIdx) + sem_dim_freq(:,snrIdx); ...
              flipud(mean_dim_freq(:,snrIdx) - sem_dim_freq(:,snrIdx))], ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = snrColors{snrIdx};
    h = [h, plot(Tlist, mean_dim_freq(:,snrIdx)', ...
                 'color', snrColors{snrIdx}, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, sprintf('snr %2.2f', snrList(snrIdx))};

end
title('Frequency domain');
xlabel('Time points per trial');
ylabel('Dimensionality');
ylim([0, 2*xDim])
legend(h, lbls, 'location', 'southoutside');
xscale('log');
axis square;
hold off;