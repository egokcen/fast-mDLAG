% gp_timescale_summary.m
%
% Description: Compare the distributions of GP timescales, across all
%              Neuropixels datasets, estimated by each method.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Choose which mDLAG-inducing results to display

NUM_INDUCE = 20;  % 20 (as in Supp. Fig. 3) or 32

load('./results/gp_timescales_mdlag_time.mat');
load('./results/gp_timescales_mdlag_freq.mat');
load(sprintf('./results/gp_timescales_smdlag_S%02d.mat', NUM_INDUCE));

%% Overlay empirical CDFs for each method

figure;
hold on;
% mDLAG-time
[f,x] = ecdf(tau_time);
h1 = plot(x, f, ...
          'color', TIMECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% mDLAG-inducing
[f,x] = ecdf(tau_sparse);
h2 = plot(x, f, ...
          'color', SPARSECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
% mDLAG-frequency
[f,x] = ecdf(tau_freq);
h3 = plot(x, f, ...
          'color', FREQCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
xlabel('GP timescale (ms)');
ylabel('Cumulative proportion of latents');
axis square;
legend([h1 h2 h3], ...
       {'Time domain', 'Inducing points', 'Freq. domain'}, ...
       'location', 'southoutside');
