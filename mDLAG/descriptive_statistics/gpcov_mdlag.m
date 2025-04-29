function k = gpcov_mdlag(params, binWidth, varargin)
%
% k = gpcov_mdlag(params, binWidth, ...)
%
% Description: Compute and (optionally) plot GP covariance functions over
%              a range of time lags specified by maxlag.
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
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
%     Optional:
%
%     showPlot -- logical; set true to plot all transformed covariance
%                 functions (default: true)
%     maxlag   -- float; maximum time lag to consider when computing
%                 covariance function; same units as binWidth.
%                 (default: 4*binWidth)
%     stepres  -- float; resolution of the computed covariance 
%                 functions, in the same units of time as binWidth
%                 (default: 0.1).
%     normalize -- logical; set true to compute cross-correlation functions,
%                  else cross-covariance functions. (default: true)
%
% Outputs:
%
%     k --(numGroups x numGroups) cell array; k{i,j} is a 
%         (1 x xDim) cell array, and k{i,j}{r} is the cross-
%         covariance function between group i and group j, given by
%         latent r. k{i,j}{r} in an anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     01 Feb 2023 -- Initial full revision.
%     02 Sep 2023 -- Overhauled to use anonymous functions, allow
%                    alternative GP kernels.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

showPlot = true;
stepres = 0.1;
maxlag = 4*binWidth;
normalize = true;
assignopts(who,varargin);

numGroups = length(params.yDims);
xDim = params.xDim;
lagSteps = -maxlag:stepres:maxlag;

% Convert GP params to same units as binWidth
gp_params = getGPparams_mdlag(params,binWidth);

% Compute prior covariance functions
k = cell(numGroups,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        k{groupIdx1,groupIdx2} = cell(1,xDim);
        if ~normalize
            CC1 = diag(sum(cat(3,params.C.moments{groupIdx1}{:}),3));
            CC2 = diag(sum(cat(3,params.C.moments{groupIdx2}{:}),3));
        end
        for xIdx = 1:xDim
            D = (gp_params.D(groupIdx2,xIdx) - gp_params.D(groupIdx1,xIdx));
            switch params.covType
                case 'rbf'
                    k{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-0.5*(gp_params.tau(xIdx)^(-2)).*(t-D).^2);
                case 'sg'
                    k{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-0.5*(gp_params.tau(xIdx)^(-2)).*(t-D).^2) ...
                          .*cos(2*pi*gp_params.nu(xIdx).*(t-D));
                case 'exp'
                    k{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-(gp_params.tau(xIdx)^(-1)).*abs(t-D));
                case 'expcos'
                    k{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-(gp_params.tau(xIdx)^(-1)).*abs(t-D)) ...
                        .*cos(2*pi*gp_params.nu(xIdx).*(t-D));
            end
            if ~normalize
                k{groupIdx1,groupIdx2}{xIdx} = @(t) sqrt(CC1(xIdx))...
                        .* k{groupIdx1,groupIdx2}{xIdx}(t) .* sqrt(CC2(xIdx));
            end
        end
    end
end

if showPlot
    for xIdx = 1:xDim
        figure;
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_curr = k{groupIdx1,groupIdx2}{xIdx}(lagSteps);
                subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                hold on;
                [maxval, maxdelay] = max(k_curr);
                minval = min(k_curr);
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [minval 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [minval 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                axis square;
                xlabel('Time lag (ms)');
                ylabel('Correlation');
                axis([min(lagSteps) max(lagSteps) minval 1.05*max([1 maxval])]);
            end
        end
    end
end