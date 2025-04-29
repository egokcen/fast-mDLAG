% bias_inducingpoints.m
%
% Description: Characterize the bias of smDLAG GP parameter estimates as a 
%              function of the number of inducing points used.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_smdlag_scaling_seqlength;
set_consts_mdlag_scaling_seqlength;

%% Characterize model performance

% Inducing points
D_sparse = nan(numS, numRuns);    % Time delays
tau_sparse = nan(numS, numRuns);  % Timescales

for runIdx = 1:numRuns
    
    fprintf('Run %d of %d...\n', runIdx, numRuns);
    
    for i = 1:numS     
        S = Slist(i);
        fprintf('    S = %d...\n', S);
        
        % Load fitted models
        results_sparse = load( ...
            sprintf('%s/run%03d/T%04d_sparse_S%03d.mat', resultDir, runIdx, Tfix, S) ...
        );
        % Estimated parameters
        % Inducing points
        gp_params_est_sparse = getGPparams_mdlag(results_sparse.estParams, binWidth);
        D_sparse(i,runIdx) = gp_params_est_sparse.D(end,1);
        tau_sparse(i,runIdx) = gp_params_est_sparse.tau(1);
        
    end

end

%% Visualize results

% Mean and SEM over runs
mean_D_sparse = mean(D_sparse,2);
sem_D_sparse = std(D_sparse,0,2) / sqrt(numRuns);
mean_tau_sparse = mean(tau_sparse,2);
sem_tau_sparse = std(tau_sparse,0,2) / sqrt(numRuns);

% Induced sampling rates
induced_rates = 1000.* Slist ./ (Tfix * binWidth);

figure;
% Timescales
subplot(1,2,1);
hold on;
% Ground truth
line([induced_rates(1) induced_rates(end)], [gp_params.tau(1) gp_params.tau(1)], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
% Inducing points
p = fill([induced_rates, fliplr(induced_rates)], ...
         [mean_tau_sparse + sem_tau_sparse; ...
          flipud(mean_tau_sparse - sem_tau_sparse)]', ...
         '', ...  % Add custom color below
         'edgecolor', 'none', ...
         'facealpha', 0.2);
p.FaceColor = SPARSECOLOR;
plot(induced_rates, mean_tau_sparse, ...
     'color', SPARSECOLOR, ...
     'linestyle', '-', ...
     'linewidth', 1.5, ...
     'marker', '.', ...
     'markersize', 10);
xlabel('Induced sampling rate (Hz)');
ylabel('Timescale (ms)');
ylim([0, gp_params.tau(1)*2.5])
axis square;
hold off;

% Time delays
subplot(1,2,2);
hold on;
% Ground truth
line([induced_rates(1) induced_rates(end)], [gp_params.D(end,1) gp_params.D(end,1)], ...
     'Color', 'k', ...
     'LineStyle', '--', ...
     'LineWidth', 1.5);
% Inducing points
p = fill([induced_rates, fliplr(induced_rates)], ...
         [mean_D_sparse + sem_D_sparse; ...
          flipud(mean_D_sparse - sem_D_sparse)]', ...
         '', ...  % Add custom color below
         'edgecolor', 'none', ...
         'facealpha', 0.2);
p.FaceColor = SPARSECOLOR;
plot(induced_rates, mean_D_sparse, ...
     'color', SPARSECOLOR, ...
     'linestyle', '-', ...
     'linewidth', 1.5, ...
     'marker', '.', ...
     'markersize', 10);
xlabel('Induced sampling rate (Hz)');
ylabel('Time delay (ms)');
ylim([9, gp_params.D(end,1)*1.1])
axis square;
hold off;