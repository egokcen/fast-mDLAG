% runtime_summary_numgroups.m
%
% Description: Quantify the runtime of mDLAG-frequency with the number of
%              groups.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Load and plot results

load('./results/runtime_mdlag_freq_numgroups.mat');
numStim = size(itertimes, 2);
numSess = size(itertimes, 3);

figure;
for sessIdx = 1:numSess
    for stimIdx = 1:numStim
        % Time per EM iteration
        subplot(1,3,1);
        hold on;
        % Reference lines
        % Linear scaling
        h2 = line([2.1 6.8], 10^(0.15).*[2.1 6.8], ...
             'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 1.5);
        % Frequency domain
        h1 = plot(groupNumbers, itertimes(:,stimIdx,sessIdx), ...
                  'color', FREQCOLOR, ...
                  'linestyle', '-', ...
                  'linewidth', 1.5, ...
                  'marker', '.', ...
                  'markersize', 10);
        xlabel('No. groups');
        ylabel('Clock time per iteration (s)');
        xscale('log');
        yscale('log');
        xlim([1.9 7.2]);
        ylim([10^(-0.5) 10^1]);
        axis square;
        legend([h1 h2], ...
               {'Freq. domain', 'Linear reference'}, ...
               'location', 'southoutside');
        hold off;
        
        % Total number of EM iterations
        subplot(1,3,2);
        hold on;
        % Frequency domain
        h1 = plot(groupNumbers, numiters(:,stimIdx,sessIdx), ...
                  'color', FREQCOLOR, ...
                  'linestyle', '-', ...
                  'linewidth', 1.5, ...
                  'marker', '.', ...
                  'markersize', 10);
        xlabel('No. groups');
        ylabel('No. iterations');
        xscale('log');
        yscale('log');
        xlim([1.9 7.2]);
        ylim([10^3 10^4.5]);
        axis square;
        legend([h1 h2], ...
               {'Freq. domain', 'Linear reference'}, ...
               'location', 'southoutside');
        hold off;
        
        % Total runtime
        subplot(1,3,3);
        hold on;
        % Reference lines
        % Linear scaling
        h2 = line([2.1 6.8], 10^(3.8).*[2.1 6.8], ...
             'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 1.5);
        % Frequency domain
        h1 = plot(groupNumbers, runtimes(:,stimIdx,sessIdx), ...
                  'color', FREQCOLOR, ...
                  'linestyle', '-', ...
                  'linewidth', 1.5, ...
                  'marker', '.', ...
                  'markersize', 10);
        xlabel('No. groups');
        ylabel('Total clock time (s)');
        xscale('log');
        yscale('log');
        xlim([1.9 7.2]);
        ylim([10^3 10^(4.8)]);
        axis square;
        legend([h1 h2], ...
               {'Freq. domain', 'Linear reference'}, ...
               'location', 'southoutside');
        hold off;

    end
end