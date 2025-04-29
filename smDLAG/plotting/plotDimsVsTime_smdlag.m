function plotDimsVsTime_smdlag(seq, xspec, params, binWidth, varargin)
%
% plotDimsVsTime_smdlag(seq, xspec, binWidth, varargin)
%
% Description: Plot each smDLAG dimension versus time in a separate panel, 
%              along with the mean trajectory across trials.
%              Group and scale panels according to which observation group
%              trajectories belong to.
%
% Arguments:
%
%     Required:
%
%     seq       -- data structure containing extracted trajectories
%     xspec     -- string; field name of trajectories in 'seq' to be 
%                  plotted (e.g., 'xorth' or 'xsm')
%     binWidth  -- bin width used when fitting model
%     params    -- Structure containing smDLAG model parameters.
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
%     Optional:
%
%     nPlotMax  -- int; maximum number of trials to plot (default: 20)
%                  NOTE: Not relevant if trialGroups is specified.
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
%     trialGroups -- (1 x numTrialGroups) cell array; Each element contains
%                    a list of trials to be grouped together for
%                    color-coding and for computing a mean time course.
%                    (default: {})
%     trialColors -- (1 x numTrialGroups) cell array; If trialGroups is
%                    specified, then trialColors must be specified, 
%                    where each element is the color for a given trial
%                    group. (default: {})
%     plotW -- logical; if true, additionally plot inducing points in a 
%              separate figure (default: false)
%     wspec -- string; field name of inducing points in 'seq' to be 
%              plotted (default: 'wsm')
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     29 Nov 2022 -- Initial full revision.

nPlotMax   = 20;
plotSingle = true;
plotMean   = true;
units      = '';
trialGroups = {};
trialColors = {};
plotW = false;
wspec = 'wsm';
assignopts(who, varargin);

colors = generateColors(); % Get custom colors for plotting
numGroups = length(params.yDims);
xDim = params.xDim;

% Check if there are any trajectories
Xall = [seq.(xspec)];

% Set up figure and axes
if isempty(Xall)
    fprintf('plotDimsVsTime_mdlag: No latents to plot.\n');
    return;
end

% Number of trajectories to be plotted
N = 0;
allTrials = [];
if isempty(trialGroups)
    N = min(length(seq), nPlotMax); 
    allTrials = 1:N;
else
    for trialGroupIdx = 1:length(trialGroups)
        N = N + length(trialGroups{trialGroupIdx}); 
        allTrials = [allTrials trialGroups{trialGroupIdx}];
    end
end
  
% Group sequences and parameters
groupSeq = partitionSeq(seq,xDim.*ones(1,numGroups),'datafield',xspec);

% Convert delay labels to units of time
gp_params = getGPparams_smdlag(params, binWidth);
D = gp_params.D;

% Global figure properties
f = figure;

% Number of rows and columns
nRows = numGroups;
nCols = xDim;

% Size and axis scales
set(f, 'Units', 'Normalized', ...
    'OuterPosition', [0.05 0.05 min([1 0.2*nCols]) min([1 0.25*nRows])]);
Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

fontsize = 12; % Size of text label fonts

for groupIdx = 1:numGroups

    % Determine vertical scale for current group
    Xall = [groupSeq{groupIdx}(allTrials).(xspec)];
    xMax = ceil(10*max(abs(Xall(:)))) / 10;
    xMid = xMax / 2; %ceil(10*(xMax/2)) / 10;
    ytk     = [-xMax -xMid 0 xMid xMax];

    for k = 1:xDim
        h = subplot(nRows, nCols, (groupIdx-1)*nCols+k);
        hold on;
        if isempty(trialGroups)
            % Initialize the mean trajectory. Only average over time points up to
            % Tmin, if trial lengths are different.
            xmean = zeros(1,Tmin); 
            for n = allTrials
                dat = groupSeq{groupIdx}(n).(xspec);
                if plotSingle 
                    % Plot single-trial trajectories
                    T = groupSeq{groupIdx}(n).T;
                    col = colors.grays{5};
                    plot(1:T, dat(k,:), ...
                         'linewidth', 0.05, ...
                         'color', col);
                end
                xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
            end
            % Plot the mean trajectory
            if plotMean
                plot(1:Tmin, xmean, ...
                     'linewidth', 2.0, ... 
                     'color', colors.grays{1});
            end     
        else
            for trialGroupIdx = 1:length(trialGroups)
                % Initialize the mean trajectory. Only average over time 
                % points up to Tmin, if trial lengths are different.
                xmean = zeros(1,Tmin); 
                for n = trialGroups{trialGroupIdx}
                    dat = groupSeq{groupIdx}(n).(xspec);
                    if plotSingle 
                        % Plot single-trial trajectories
                        T = groupSeq{groupIdx}(n).T;
                        plot(1:T, dat(k,:), ...
                             'linewidth', 0.05, ...
                             'color', trialColors{trialGroupIdx});
                    end
                    xmean = xmean + (1.0/length(trialGroups{trialGroupIdx})) ...
                            .* dat(k,1:Tmin);
                end
                % Plot the mean trajectory
                if plotMean
                    plot(1:Tmin, xmean, ...
                         'linewidth', 2.0, ... 
                         'color', trialColors{trialGroupIdx});
                end
                
            end
        end
        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}_{%d,%d,:}$$',groupIdx,k);
        str = sprintf('%s, $$D_{%d,%d} = %3.1f$$', str, groupIdx, k, D(groupIdx,k));
        str = sprintf('%s%s', str, units);
        title(str, 'interpreter', 'latex', 'fontsize', fontsize);
        xlabel(sprintf('Time%s', units));
    end

end

if plotW
    
    f = figure;

    % Size and axis scales
    set(f, 'Units', 'Normalized', ...
        'OuterPosition', [0.05 0.05 min([1 0.2*nCols]) min([1 0.25])]);
    
    % Determine vertical scale
    Wall = [];
    for n = allTrials
        for k = 1:xDim
            Wall = [Wall seq(n).(wspec){k}]; 
        end
    end
    wMax = ceil(10*max(abs(Wall))) / 10;
    wMid = wMax / 2;
    ytk  = [-wMax -wMid 0 wMid wMax];

    for k = 1:xDim
        Sk = length(params.Z{k});
        [Zk,sortIdxs] = sort(params.Z{k});
        h = subplot(1, nCols, k);
        hold on;
        if isempty(trialGroups)
            % Initialize the mean trajectory. Only average over time points up to
            % Tmin, if trial lengths are different.
            wmean = zeros(1,Sk); 
            for n = allTrials
                if plotSingle 
                    % Plot single-trial trajectories
                    col = colors.grays{5};
                    plot(Zk, seq(n).(wspec){k}(sortIdxs), ...
                         'linewidth', 0.05, ...
                         'color', col);
                end
                wmean = wmean + (1.0/N) .* seq(n).(wspec){k};
            end
            % Plot the mean trajectory
            if plotMean
                plot(Zk, wmean(sortIdxs), ...
                     'linewidth', 2.0, ... 
                     'color', colors.grays{1});
            end     
        else
            for trialGroupIdx = 1:length(trialGroups)
                % Initialize the mean trajectory. Only average over time 
                % points up to Tmin, if trial lengths are different.
                wmean = zeros(1,Sk); 
                for n = trialGroups{trialGroupIdx}
                    if plotSingle 
                        % Plot single-trial trajectories
                        plot(Zk, seq(n).(wspec){k}(sortIdxs), ...
                             'linewidth', 0.05, ...
                             'color', trialColors{trialGroupIdx});
                    end
                    xmean = xmean + (1.0/length(trialGroups{trialGroupIdx})) ...
                            .* seq(n).(wspec){k};
                end
                % Plot the mean trajectory
                if plotMean
                    plot(Zk, wmean(sortIdxs), ...
                         'linewidth', 2.0, ... 
                         'color', trialColors{trialGroupIdx});
                end

            end
        end
        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf w}_{%d,:}$$',k);
        title(str, 'interpreter', 'latex', 'fontsize', fontsize);
        xlabel(sprintf('Time%s', units));
    end 
    
end