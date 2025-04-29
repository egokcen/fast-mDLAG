function plotFittingProgress(trackedParams, binWidth, covType, varargin)
%
% plotFittingProgress(trackedParams, binWidth, covType, varargin)
%
% Description: Plot intermediate values stored throughout mDLAG model
%              fitting vs iteration number. Visualizing these values can
%              help with troubleshooting.
%
% Arguments:
%
%     Required:
%
%     trackedParams -- structure containing parameters tracked  
%                      throughout fitting:
%         lb       -- (1 x numIters) array; variational lower bound 
%                      at each iteration.
%                      Note: Entries will be NaN, on iterations
%                      was lb was not computed, to save time.
%         iterTime -- (1 x numIters) array; computation time for 
%                     each EM iteration
%         alphas   -- (1 x numIters) cell arry; estimated ARD 
%                     parameters (alpha.mean) after each EM
%                     iteration.
%         gp_params -- structure tracking kernel-dependent GP parameters. 
%                      Contains the following fields:
%           For 'rbf' or 'exp':
%             Ds       -- (1 x numIters) cell array; the estimated delay 
%                         matrix (D) after each EM iteration.
%             gams     -- (1 x numIters) cell arry; estimated gamma after 
%                         each EM iteration.
%           For 'sg' or 'expcos':
%             Ds       -- (1 x numIters) cell array; the estimated delay 
%                         matrix (D) after each EM iteration.
%             gams     -- (1 x numIters) cell arry; estimated gamma after 
%                         each EM iteration.
%             nus      -- (1 x numIters) cell arry; estimated nu after
%                         each EM iteration.
%     binWidth -- float; bin width or sample period, in units of time 
%                 (e.g., ms)
%     covType  -- string; type of GP covariance. The following GP kernels 
%                 are currently supported:
%                     'rbf'    -- Radial basis function, or squared 
%                                 exponential
%                     'sg'     -- Spectral Gaussian
%                     'exp'    -- Exponential
%                     'expcos' -- Exponential-cosine
%
%     Optional:
%     
%     freqLB    -- int; if known, specify how often LL was computed during
%                  model fitting (default: 1)
%     freqParam -- int; if known, specify how often GP parameters were
%                  stored during model fitting (default: 1)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
%
% Outputs:
%     None. (But creates a figure)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Oct 2022 -- Initial full revision.
%     08 Dec 2022 -- Changed input argument structure.
%     02 Sep 2023 -- Added ability to handle alternate GP kernels.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

freqLB = 1;
freqParam = 1;
units = '';
assignopts(who,varargin);

% Format units for axis labels
if ~isempty(units)
    timeUnits = sprintf(' (%s)', units); 
    freqUnits = sprintf(' (1/%s)', units);
end

[numGroups,xDim] = size(trackedParams.gp_params.Ds{1});

colors = generateColors(); % Get custom colors for plotting
fontsize = 12;

figure;

% Lower bound curve
subplot(1,2,1);
hold on;
lb = trackedParams.lb(~isnan(trackedParams.lb));
xlb = [];
if freqLB <= 1
    xlb = 1:length(lb);
elseif freqLB == 2
    xlb = [1 freqLB.*(1:length(lb)-1)];
else
    xlb = [1 2 freqLB.*(1:length(lb)-2)];
end
plot(xlb, lb, 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Lower bound');

% Cumulative fitting time
subplot(1,2,2);
hold on;
plot(cumsum(trackedParams.iterTime), 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Cumulative time elapsed (s)');

if xDim > 0
    if numGroups > 1
        % Delay progress
        figure;
        Delays = cat(3, trackedParams.gp_params.Ds{:});
        for rowIdx = 1:numGroups
            for colIdx = 1:numGroups
                % Avoid redundant plots, and plot only different groups
                if colIdx > rowIdx
                    subplot(numGroups,numGroups,numGroups*(rowIdx-1)+colIdx);
                    hold on;
                    for i = 1:xDim
                        plot(freqParam.*(0:length(trackedParams.gp_params.Ds)-1), ...
                             binWidth.*squeeze(Delays(colIdx,i,:) - Delays(rowIdx,i,:)),...
                             'linewidth', 1.5);
                    end
                    xlabel('# Iterations');
                    ylabel(sprintf('Delay%s, group %d to %d', timeUnits, rowIdx, colIdx));
                end
            end
        end
    end

    if ismember(covType, {'rbf', 'sg'})
        % GP timescale progress
        figure; 
        hold on;
        Gams = cat(3, trackedParams.gp_params.gams{:});
        for i = 1:xDim
            plot(freqParam.*(0:length(trackedParams.gp_params.gams)-1), ...
                 binWidth./sqrt(squeeze(Gams(1,i,:))), ...
                 'linewidth', 1.5);
        end
        xlabel('# Iterations');
        ylabel(sprintf('GP timescale%s', timeUnits));
    end

    if ismember(covType, {'exp', 'expcos'})
        % GP timescale progress
        figure; 
        hold on;
        Gams = cat(3, trackedParams.gp_params.gams{:});
        for i = 1:xDim
            plot(freqParam.*(0:length(trackedParams.gp_params.gams)-1), ...
                 binWidth./squeeze(Gams(1,i,:)), ...
                 'linewidth', 1.5);
        end
        xlabel('# Iterations');
        ylabel(sprintf('GP timescale%s', timeUnits));
    end

    if ismember(covType, {'sg', 'expcos'})
        % GP center frequency progress
        figure; 
        hold on;
        Nus = cat(3, trackedParams.gp_params.nus{:});
        for i = 1:xDim
            plot(freqParam.*(0:length(trackedParams.gp_params.nus)-1), ...
                 squeeze(Nus(1,i,:))./binWidth, ...
                 'linewidth', 1.5);
        end
        xlabel('# Iterations');
        ylabel(sprintf('GP frequencies%s', freqUnits));
    end
    
end

% ARD parameter progress
alphas = cat(3, trackedParams.alphas{:});
figure;
for groupIdx = 1:numGroups
    subplot(1,numGroups,groupIdx);
    hold on;
    for i = 1:xDim
        plot(freqParam.*(0:length(trackedParams.alphas)-1), ...
             1./squeeze(alphas(groupIdx,i,:)), ...
             'linewidth', 1.5);
    end
    xlabel('# Iterations');
    ylabel('\alpha^{-1}', 'fontsize', fontsize);
    title(sprintf('Group %d', groupIdx));
end