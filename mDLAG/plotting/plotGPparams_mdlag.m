function gp_params = plotGPparams_mdlag(params,binWidth,varargin)
%
% gp_params = plotGPparams_mdlag(params,binWidth,...)
%
% Description: Plot mDLAG Delay Matrix and latent timescales. Delays and 
%              timescales are plotted as ordered pairs.
%
% Arguments: 
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
%         covType -- string; type of GP covariance. The 
%                    following GP kernels are currently 
%                    supported:
%                        'rbf'    -- Radial basis function, or 
%                                    squared exponential
%                        'sg'     -- Spectral Gaussian
%                        'exp'    -- Exponential
%                        'expcos' -- Exponential-cosine
%         The following fields depend on GP covariance type:
%             For 'rbf':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)                                                      
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'sg':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)    
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'exp':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma                                                  
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'expcos':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma    
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
%         C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                       mean loadings matrix for each group
%         C.covs     -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       covariance of a row of C.
%         C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       second moment of a row of C.
%         alpha.a    -- (numGroups x 1) array; shape parameters of alpha 
%                       posterior
%         alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                       alpha posterior
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights (for ARD); alpha.a ./ alpha.b
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%
%     Optional:
%
%     sigDims  -- (numGroups x xDim) logical array; sigDims(i,j) is 1 if
%                 latent j explains significant shared variance in group i,
%                 0 otherwise. (default: [])
%     units    -- string; units of time of binWidth (for labels)
%                 (default: '')
%
% Outputs:
%     
%    gp_params -- structure containing mDLAG GP parameters, converted into
%                 units of time (or 1/time).
%                 For 'rbf' or 'exp':
%                     D   -- (numGroups x xDim) array; delays from latents
%                            to observed variables
%                     tau -- (1 x xDim) array; GP timescales  
%                 For 'sg' or 'expcos':
%                     D   -- (numGroups x xDim) array; delays from latents
%                            to observed variables
%                     tau -- (1 x xDim) array; GP timescales  
%                     nu  -- (1 x xDim) array; center frequencies
%              
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Oct 2022 -- Initial full revision.
%     01 Sep 2023 -- Overhauled for better GP kernel modularity.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

% Set optional arguments
sigDims = [];
units = '';
assignopts(who,varargin);

xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
if isempty(sigDims)
    % Plot GP parameters for all latents in all groups, regardless of 
    % significance
    sigDims = ones(numGroups,xDim);
end

colors = generateColors(); % Generate custom plotting colors
pointsize = 25; % Size of scatterplot points
% Format units for axis labels
if ~isempty(units)
    timeUnits = sprintf(' (%s)', units); 
    freqUnits = sprintf(' (1/%s)', units);
end
nyqFreq = 0.5/binWidth; % Maximum samplable frequency without aliasing

% Convert GP params into units of time
gp_params = getGPparams_mdlag(params, binWidth);
maxTau = max(gp_params.tau);

% Plot shared latents and pairwise delays.
% But don't make a plot if there are no shared latents.
plotShared = false;
for i = 1:numGroups
    for j = 1:numGroups
        % Consider only different groups
        if i ~= j
            % Determine which latents are shared between the current pair
            % of groups
            sig = find(all(sigDims([i j],:),1));
            if ~isempty(sig)
                plotShared = true; 
            end
        end
    end
end

if plotShared
    figure;
    for i = 1:numGroups
        for j = 1:numGroups

            % Consider only different groups
            if i ~= j

                % Take only the delays associated with the current pair of groups
                delays = gp_params.D(j,:) ...
                       - gp_params.D(i,:);

                % Determine which latents are shared between the current pair of
                % groups
                sig = find(all(sigDims([i j],:),1));

                % Figure out horizontal axis limits
                maxDelay = max(abs(delays));
                if maxDelay <= 0.05*maxTau
                    maxDelay = 0.5*maxTau;
                end

                switch params.covType

                    case {'rbf', 'exp'}
                        
                        subplot(numGroups,numGroups,(i-1)*numGroups+j);
                        hold on;
                        xlabel(sprintf('Delay from group %d to group %d%s',i,j,timeUnits));
                        ylabel(sprintf('GP timescale%s', timeUnits));
                        xlim([-1.1*maxDelay,1.1*maxDelay]);
                        ylim([0,1.1*maxTau]);
                        line([0 0], [0 1.1*maxTau], ...
                             'Color', colors.grays{6}, ...
                             'linestyle', '--', ...
                             'linewidth', 1.5);
                        plot(delays(sig), gp_params.tau(sig), ...
                             'marker', '.', ...
                             'linestyle', 'none', ...
                             'color', colors.grays{1}, ...
                             'markersize', pointsize);
                        hold off;

                    case {'sg', 'expcos'}
        
                        % Figure out vertical axis limits
                        maxY = max(gp_params.nu + 2./(2*pi*gp_params.tau));
                        maxY = max([maxY nyqFreq]);
                        minY = min(gp_params.nu - 2./(2*pi*gp_params.tau));
                        subplot(numGroups,numGroups,(i-1)*numGroups+j);
                        hold on;
                        xlabel(sprintf('Delay from group %d to group %d%s',i,j,timeUnits));
                        ylabel(sprintf('GP frequency%s', freqUnits));
                        xlim([-1.1*maxDelay,1.1*maxDelay]);
                        ylim([minY,maxY]);
                        line([0 0], [minY maxY], ...
                             'Color', colors.grays{6}, ...
                             'linestyle', '--', ...
                             'linewidth', 1.5);
                        line([-1.1*maxDelay,1.1*maxDelay], [nyqFreq nyqFreq], ...
                             'Color', colors.reds{6}, ...
                             'linestyle', '--', ...
                             'linewidth', 1.5);
                        errorbar(delays(sig), gp_params.nu(sig), 2./(2*pi*gp_params.tau(sig)), ...
                                 'marker', 'none', ...
                                 'linestyle', 'none', ...
                                 'linewidth', 1.5, ...
                                 'color', colors.grays{1});
                        plot(delays(sig), gp_params.nu(sig), ...
                             'marker', '.', ...
                             'linestyle', 'none', ...
                             'color', colors.grays{1}, ...
                             'markersize', pointsize);
                        hold off;
                end
            end
        end
    end
end


% Plot local latents (i.e., latents present in only one group)
% First determine how many local latents are in each group, and
% don't attempt to plot if there are no local latents
localDims = cell(1,numGroups);
plotLocal = false;
for groupIdx = 1:numGroups
    local = zeros(numGroups,1);
    local(groupIdx) = 1;
    localDims{groupIdx} = find(ismember(sigDims',local','rows'));
    if ~isempty(localDims{groupIdx})
        plotLocal = true;
    end
end

if plotLocal
    figure;
    for groupIdx = 1:numGroups
        numLocal = length(localDims{groupIdx});
        % Don't try to plot anything for groups with 0 local latents
        if numLocal > 0
            switch params.covType

                case {'rbf', 'exp'}

                    subplot(1,numGroups,groupIdx);
                    hold on;
                    xlim([0,numLocal+1]); 
                    ylim([0,1.1*maxTau]);
                    h = bar(1:numLocal,gp_params.tau(localDims{groupIdx}),0.4);    
                    set(h,'facecolor',colors.grays{1},'edgecolor',colors.grays{1});
                    ylabel(sprintf('GP timescale%s', timeUnits));
                    set(gca,'XTick',1:numLocal);
                    xlabel(sprintf('Local latents, group %d', groupIdx)); 
                    hold off;

                case {'sg', 'expcos'}

                    maxY = max(gp_params.nu(localDims{groupIdx}) ...
                        + 2./(2*pi*gp_params.tau(localDims{groupIdx})));
                    maxY = max([maxY nyqFreq]);
                    minY = min(gp_params.nu(localDims{groupIdx}) ...
                        - 2./(2*pi*gp_params.tau(localDims{groupIdx})));
                    subplot(1,numGroups,groupIdx);
                    hold on;
                    xlim([0,numLocal+1]);
                    ylim([minY maxY]);
                    line([-10,+10], [nyqFreq nyqFreq], ...
                         'Color', colors.reds{6}, ...
                         'linestyle', '--', ...
                         'linewidth', 1.5);
                    errorbar(1:numLocal, gp_params.nu(localDims{groupIdx}), ...
                             2./(2*pi*gp_params.tau(localDims{groupIdx})), ...
                             'marker', 'none', ...
                             'linestyle', 'none', ...
                             'linewidth', 1.5, ...
                             'color', colors.grays{1});
                    plot(1:numLocal,gp_params.nu(localDims{groupIdx}), ...
                         'marker', '.', ...
                         'linestyle', 'none', ...
                         'color', colors.grays{1}, ...
                         'markersize', pointsize);
                    ylabel(sprintf('GP frequency%s', freqUnits));
                    set(gca,'XTick',1:numLocal);
                    xlabel(sprintf('Local latents, group %d', groupIdx)); 
                    hold off;
            end
        end
    end
end
