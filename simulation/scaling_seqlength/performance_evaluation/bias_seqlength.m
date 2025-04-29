% bias_seqlength.m
%
% Description: Characterize the bias of mDLAG GP parameter estimates as a 
%              function of sequence (trial) length. Compare time, inducing, 
%              and frequency domain fitting approaches. Optionally display
%              results from extra experiments, including tapering and
%              fixing certain GP parameters
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_scaling_seqlength;

%% Characterize model performance

FIX_TYPE = 0;  % 0: no parameter fixing; 1: fix timescales; 2: fix delays
PLOT_SPARSE = 0; % 0: skip smDLAG results; 1: plot smDLAG results
PLOT_TAPER = 1;  % 0: skip taper results; 1: plot taper results

if ~FIX_TYPE
    % Time domain
    D_time = nan(numPartitions, numRuns);    % Time delays
    tau_time = nan(numPartitions, numRuns);  % Timescales

    if PLOT_SPARSE
        % Inducing points
        D_sparse = nan(numPartitions, numRuns);    % Time delays
        tau_sparse = nan(numPartitions, numRuns);  % Timescales
    end

    if PLOT_TAPER
        % Frequency domain with taper
        D_taper = nan(numPartitions, numRuns);     % Time delays
        tau_taper = nan(numPartitions, numRuns);   % Timescales
    end
end

% Frequency domain
D_freq = nan(numPartitions, numRuns);    % Time delays
tau_freq = nan(numPartitions, numRuns);  % Timescales

for runIdx = 1:numRuns
    
    fprintf('Run %d of %d...\n', runIdx, numRuns);
    
    for i = 1:numPartitions
        T = Tlist(i);              
        fprintf('    T = %d...\n', T);
        
        % Load fitted models
        switch FIX_TYPE
            case 0
                results_time = load( ...
                    sprintf('%s/run%03d/T%04d_time.mat', resultDir, runIdx, T) ...
                );
                if PLOT_SPARSE
                    results_sparse = load( ...
                        sprintf('%s/run%03d/T%04d_sparse.mat', resultDir, runIdx, T) ...
                    );
                end
                results_freq = load( ...
                    sprintf('%s/run%03d/T%04d_freq.mat', resultDir, runIdx, T) ...
                );
                if PLOT_TAPER
                    results_taper = load( ...
                        sprintf('%s/run%03d/T%04d_freq_hamming.mat', resultDir, runIdx, T) ...
                    );
                end
            case 1
                results_freq = load( ...
                    sprintf('%s/run%03d/T%04d_freq_fixtau.mat', resultDir, runIdx, T) ...
                );
            case 2
                results_freq = load( ...
                    sprintf('%s/run%03d/T%04d_freq_fixD.mat', resultDir, runIdx, T) ...
                );
        end
        % Estimated parameters
        if ~FIX_TYPE
            % Time domain
            gp_params_est_time = getGPparams_mdlag(results_time.estParams, binWidth);
            D_time(i,runIdx) = gp_params_est_time.D(end,1);
            tau_time(i,runIdx) = gp_params_est_time.tau(1);
            if PLOT_SPARSE
                % Inducing points
                gp_params_est_sparse = getGPparams_mdlag(results_sparse.estParams, binWidth);
                D_sparse(i,runIdx) = gp_params_est_sparse.D(end,1);
                tau_sparse(i,runIdx) = gp_params_est_sparse.tau(1);
            end
            if PLOT_TAPER
                % Frequency domain with taper
                gp_params_est_taper = getGPparams_mdlag(results_taper.estParams, binWidth);
                D_taper(i,runIdx) = gp_params_est_taper.D(end,1);
                tau_taper(i,runIdx) = gp_params_est_taper.tau(1);
            end
        end

        % Frequency domain
        gp_params_est_freq = getGPparams_mdlag(results_freq.estParams, binWidth);
        D_freq(i,runIdx) = gp_params_est_freq.D(end,1);
        tau_freq(i,runIdx) = gp_params_est_freq.tau(1);
        
    end

end

%% Visualize results

% Mean and SEM over runs
if ~FIX_TYPE
    % Time domain
    mean_D_time = mean(D_time,2);
    sem_D_time = std(D_time,0,2) / sqrt(numRuns);
    mean_tau_time = mean(tau_time,2);
    sem_tau_time = std(tau_time,0,2) / sqrt(numRuns);

    if PLOT_SPARSE
        % Inducing points
        mean_D_sparse = mean(D_sparse,2);
        sem_D_sparse = std(D_sparse,0,2) / sqrt(numRuns);
        mean_tau_sparse = mean(tau_sparse,2);
        sem_tau_sparse = std(tau_sparse,0,2) / sqrt(numRuns);
    end

    if PLOT_TAPER
        % Frequency domain with taper
        mean_D_taper = mean(D_taper,2);
        sem_D_taper = std(D_taper,0,2) / sqrt(numRuns);
        mean_tau_taper = mean(tau_taper,2);
        sem_tau_taper = std(tau_taper,0,2) / sqrt(numRuns);
    end
end

% Frequency domain
mean_D_freq = mean(D_freq,2);
sem_D_freq = std(D_freq,0,2) / sqrt(numRuns);
mean_tau_freq = mean(tau_freq,2);
sem_tau_freq = std(tau_freq,0,2) / sqrt(numRuns);

figure;
% Timescales
h = [];
lbls = {};
subplot(1,2,1);
hold on;
% Ground truth
line([Tlist(1) Tlist(end)], [gp_params.tau(1) gp_params.tau(1)], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
if ~FIX_TYPE
    % Time domain
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_tau_time + sem_tau_time; ...
              flipud(mean_tau_time - sem_tau_time)]', ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = TIMECOLOR;
    h = [h, plot(Tlist, mean_tau_time, ...
                 'color', TIMECOLOR, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, 'Time domain'};
    if PLOT_SPARSE
        % Inducing points
        p = fill([Tlist, fliplr(Tlist)], ...
                 [mean_tau_sparse + sem_tau_sparse; ...
                  flipud(mean_tau_sparse - sem_tau_sparse)]', ...
                 '', ...  % Add custom color below
                 'edgecolor', 'none', ...
                 'facealpha', 0.2);
        p.FaceColor = SPARSECOLOR;
        h = [h, plot(Tlist, mean_tau_sparse, ...
                     'color', SPARSECOLOR, ...
                     'linestyle', '-', ...
                     'linewidth', 1.5, ...
                     'marker', '.', ...
                     'MarkerSize', 10)];
        lbls = {lbls{:}, 'Inducing points'};
    end
end
% Frequency domain
p = fill([Tlist, fliplr(Tlist)], ...
         [mean_tau_freq + sem_tau_freq; ...
          flipud(mean_tau_freq - sem_tau_freq)]', ...
         '', ...  % Add custom color below
         'edgecolor', 'none', ...
         'facealpha', 0.2);
p.FaceColor = FREQCOLOR;
h = [h, plot(Tlist, mean_tau_freq, ...
             'color', FREQCOLOR, ...
             'linestyle', '-', ...
             'linewidth', 1.5, ...
             'marker', '.', ...
             'markersize', 10)];
lbls = {lbls{:}, 'Freq. domain'};
if PLOT_TAPER && ~FIX_TYPE
    % Frequency domain with taper
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_tau_taper + sem_tau_taper; ...
              flipud(mean_tau_taper - sem_tau_taper)]', ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = TAPERCOLOR;
    h = [h, plot(Tlist, mean_tau_taper, ...
                 'color', TAPERCOLOR, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, 'Taper'};
end
xlabel('Time points per trial');
ylabel('Timescale (ms)');
ylim([0, gp_params.tau(1)*1.2])
axis square;
legend(h, lbls, 'location', 'southoutside');
hold off;

% Time delays
h = [];
lbls = {};
subplot(1,2,2);
hold on;
% Ground truth
line([Tlist(1) Tlist(end)], [gp_params.D(end,1) gp_params.D(end,1)], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
if ~FIX_TYPE
    % Time domain
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_D_time + sem_D_time; ...
              flipud(mean_D_time - sem_D_time)]', ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = TIMECOLOR;
    h = [h, plot(Tlist, mean_D_time, ...
                 'color', TIMECOLOR, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'MarkerSize', 10)];
    lbls = {lbls{:}, 'Time domain'};
    if PLOT_SPARSE
        % Inducing points
        p = fill([Tlist, fliplr(Tlist)], ...
                 [mean_D_sparse + sem_D_sparse; ...
                  flipud(mean_D_sparse - sem_D_sparse)]', ...
                 '', ...  % Add custom color below
                 'edgecolor', 'none', ...
                 'facealpha', 0.2);
        p.FaceColor = SPARSECOLOR;
        h = [h, plot(Tlist, mean_D_sparse, ...
                     'color', SPARSECOLOR, ...
                     'linestyle', '-', ...
                     'linewidth', 1.5, ...
                     'marker', '.', ...
                     'MarkerSize', 10)];
        lbls = {lbls{:}, 'Inducing points'};
    end
end
% Frequency domain
p = fill([Tlist, fliplr(Tlist)], ...
         [mean_D_freq + sem_D_freq; ...
          flipud(mean_D_freq - sem_D_freq)]', ...
         '', ...  % Add custom color below
         'edgecolor', 'none', ...
         'facealpha', 0.2);
p.FaceColor = FREQCOLOR;
h = [h, plot(Tlist, mean_D_freq, ...
             'color', FREQCOLOR, ...
             'linestyle', '-', ...
             'linewidth', 1.5, ...
             'marker', '.', ...
             'markersize', 10)];
lbls = {lbls{:}, 'Freq. domain'};
if PLOT_TAPER && ~FIX_TYPE
    % Frequency domain with taper
    p = fill([Tlist, fliplr(Tlist)], ...
             [mean_D_taper + sem_D_taper; ...
              flipud(mean_D_taper - sem_D_taper)]', ...
             '', ...  % Add custom color below
             'edgecolor', 'none', ...
             'facealpha', 0.2);
    p.FaceColor = TAPERCOLOR;
    h = [h, plot(Tlist, mean_D_taper, ...
                 'color', TAPERCOLOR, ...
                 'linestyle', '-', ...
                 'linewidth', 1.5, ...
                 'marker', '.', ...
                 'markersize', 10)];
    lbls = {lbls{:}, 'Taper'};
end
xlabel('Time points per trial');
ylabel('Time delay (ms)');
ylim([0, gp_params.D(end,1)*1.2]);
axis square;
legend(h, lbls, 'location', 'southoutside');
hold off;