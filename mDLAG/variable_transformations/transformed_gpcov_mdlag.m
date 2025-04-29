function k = transformed_gpcov_mdlag(V, params, binWidth, varargin)
%
% k = transformed_gpcov_mdlag(V, params, binWidth, ...)
%
% Description: Given a coupled set of basis vectors, V, compute GP
%              covariance functions for projections that lie along those
%              directions.
%
% Arguments:
%
%     Required:
%
%     V       -- (1 x numGroups) cell array; V{i} is a (yDims(i) x r)
%                array, where r is the number of basis vectors.
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
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
%     Optional:
%
%     showPlot  -- logical; set true to plot all transformed covariance
%                  functions (default: true)
%     maxlag    -- float; maximum time lag to consider when computing
%                  covariance function; same units as binWidth.
%                  (default: 4*binWidth)
%     stepres   -- float; resolution of the computed covariance 
%                  functions, in the same units of time as binWidth
%                  (default: 0.1).
%     normalize -- logical; true to compute correlation functions (i.e.,
%                  normalized); false to compute covariance functions
%                  (default: true)
%
% Outputs:
%
%     k -- (1 x r) cell array; k{r} is a (numGroups x numGroups) cell
%          array, and k{r}{i,j} is the cross-covariance function between
%          group i and group j. k{r}{i,j} is an anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

showPlot = true;
stepres = 0.1;
maxlag = 4*binWidth;
normalize = true;
assignopts(who,varargin);

numGroups = length(V);
r = size(V{1},2);
xDim = params.xDim;
lagSteps = -maxlag:stepres:maxlag;

% Initialize output structure
k = cell(1,r);
for j = 1:r
    k{j} = cell(numGroups); 
end

% Compute prior covariance functions
k_prior = gpcov_mdlag(params, binWidth, ...
                      'normalize', true, ...
                      'showPlot', false);

% Compute transformed covariance functions
for j = 1:r
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            k{j}{groupIdx1,groupIdx2} = @(t) 0;
            for xIdx = 1:xDim
                k{j}{groupIdx1,groupIdx2} = @(t) k{j}{groupIdx1,groupIdx2}(t) ...
                    + k_prior{groupIdx1,groupIdx2}{xIdx}(t) ...
                    .*(V{groupIdx1}(:,j)'*params.C.means{groupIdx1}(:,xIdx)...
                    *params.C.means{groupIdx2}(:,xIdx)'*V{groupIdx2}(:,j));
            end
        end
    end
end

if normalize
    % Normalize covariance functions to get correlation functions
    k_norm = k;
    for j = 1:r
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_norm{j}{groupIdx1,groupIdx2} ...
                    = @(t) k{j}{groupIdx1,groupIdx1}(0).^(-0.5) ...
                      .* k{j}{groupIdx1,groupIdx2}(t) ...
                      .* k{j}{groupIdx2,groupIdx2}(0).^(-0.5);
            end
        end
    end
    k = k_norm;
end

if showPlot
    if normalize
        ylbl = 'Correlation';
    else
        ylbl = 'Covariance';
    end
    
    for j = 1:r
        figure;
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_curr = k{j}{groupIdx1,groupIdx2}(lagSteps);
                subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                hold on;
                [maxval, maxdelay] = max(abs(k_curr));
                if k_curr(maxdelay) < 0
                    k_curr = -1.*k_curr; 
                end
                minval = min(k_curr);
                minval = 0;
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [minval 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [minval 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                axis square;
                xlabel('Time lag (ms)');
                ylabel(ylbl);
                axis([min(lagSteps) max(lagSteps) minval 1.05*max([1 maxval])]);
            end
        end
    end
end