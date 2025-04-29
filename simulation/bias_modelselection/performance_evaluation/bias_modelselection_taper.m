% bias_modelselection_taper.m
%
% Description: Characterize the bias of mDLAG dimensionality estimates as a 
%              function of sequence (trial) length. Investigate the effects
%              of tapering schemes on that bias.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_modelselection;

%% Characterize estimates of dimensionality
cutoff_sharedvar = 0.01; % Minimum shared variance within a group that a latent must explain
cutoff_snr = 0.001;      % Minimum SNR that a group must have for latents to be signficiant
SIG_IDX = 2;             % Access significant dimensions for one group

% No taper
dim_notaper = nan(numPartitions, snrPartitions, runsPerSNR);

% Hamming
dim_hamming = nan(numPartitions, snrPartitions, runsPerSNR);

for snrIdx = 1:snrPartitions
 
    for runIdx = 1:runsPerSNR

        fprintf('Run %d of %d\n', runIdx+(snrIdx-1)*runsPerSNR, numRuns);
    
        for i = 1:numPartitions
            T = Tlist(i);

            % Load fitted models
            results_notaper = load( ...
                sprintf('%s/run%03d/T%04d_freq.mat', ...
                        resultDir, ...
                        runIdx+(snrIdx-1)*runsPerSNR, ...
                        T) ...
            );
            results_hamming = load( ...
                sprintf('%s/run%03d/T%04d_freq_hamming.mat', ...
                        resultDir, ...
                        runIdx+(snrIdx-1)*runsPerSNR, ...
                        T) ...
            );

            % Estimated dimensionalities
            % No taper
            [dims,~,~,~] ...
                = computeDimensionalities(results_notaper.estParams, ...
                                          cutoff_sharedvar, ...
                                          cutoff_snr);
            dim_notaper(i,snrIdx,runIdx) = dims(SIG_IDX);
            % Hamming
            [dims,~,~,~] ...
                = computeDimensionalities(results_hamming.estParams, ...
                                          cutoff_sharedvar, ...
                                          cutoff_snr);
            dim_hamming(i,snrIdx,runIdx) = dims(SIG_IDX);

        end

    end

end

%% Visualize results

% Mean and SEM over runs
% No taper
mean_dim_notaper = mean(dim_notaper,3);
sem_dim_notaper = std(dim_notaper,0,3) / sqrt(runsPerSNR);

% Hamming
mean_dim_hamming = mean(dim_hamming,3);
sem_dim_hamming = std(dim_hamming,0,3) / sqrt(runsPerSNR);

% No taper
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
             [mean_dim_notaper(:,snrIdx) + sem_dim_notaper(:,snrIdx); ...
              flipud(mean_dim_notaper(:,snrIdx) - sem_dim_notaper(:,snrIdx))], ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = snrColors{snrIdx};
    h = [h, plot(Tlist, mean_dim_notaper(:,snrIdx)', ...
                 'color', snrColors{snrIdx}, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, sprintf('snr %2.2f', snrList(snrIdx))};

end
title('No taper');
xlabel('Time points per trial');
ylabel('Dimensionality');
ylim([0, 2*xDim])
legend(h, lbls, 'location', 'southoutside');
xscale('log');
axis square;
hold off;

% Hamming
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
             [mean_dim_hamming(:,snrIdx) + sem_dim_hamming(:,snrIdx); ...
              flipud(mean_dim_hamming(:,snrIdx) - sem_dim_hamming(:,snrIdx))], ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = snrColors{snrIdx};
    h = [h, plot(Tlist, mean_dim_hamming(:,snrIdx)', ...
                 'color', snrColors{snrIdx}, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, sprintf('snr %2.2f', snrList(snrIdx))};

end
title('Hamming');
xlabel('Time points per trial');
ylabel('Dimensionality');
ylim([0, 2*xDim])
legend(h, lbls, 'location', 'southoutside');
xscale('log');
axis square;
hold off;