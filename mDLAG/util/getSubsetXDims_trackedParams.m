function [gp_params,alphas] = getSubsetXDims_trackedParams(covType,gp_params,alphas,xDims)
%
% [gp_params,alphas] = getSubsetXDims_trackedParams(covType,gp_params,alphas,xDims)
%
% Description: Remove unwanted latents from mDLAG model parameters that
%              are tracked throughout fitting.
%
% Arguments:
%
%     covType  -- string; type of GP covariance. The following GP kernels 
%                 are currently supported:
%                     'rbf'    -- Radial basis function, or squared 
%                                 exponential
%                     'sg'     -- Spectral Gaussian
%                     'exp'    -- Exponential
%                     'expcos' -- Exponential-cosine
%     gp_params -- structure tracking kernel-dependent GP parameters. 
%                  Contains the following fields:
%         For 'rbf' or 'exp':
%           Ds       -- (1 x numIters) cell array; the estimated delay 
%                       matrix (D) after each EM iteration.
%           gams     -- (1 x numIters) cell arry; estimated gamma after 
%                       each EM iteration.
%         For 'sg' or 'expcos':
%           Ds       -- (1 x numIters) cell array; the estimated delay 
%                       matrix (D) after each EM iteration.
%           gams     -- (1 x numIters) cell arry; estimated gamma after 
%                       each EM iteration.
%           nus      -- (1 x numIters) cell arry; estimated nu after
%                       each EM iteration.
%     alphas    -- (1 x numIters) cell arry; estimated ARD parameters
%                  (alpha.mean) after each EM iteration.
%     xDims     -- (1 x numDims) array; latent state dimensions to be  
%                  retained in outparams.
%
% Outputs:
%
%     gp_params -- same structure as gp_params above, but with relevant
%                  xDims removed.
%     alphas    -- (1 x numIters) cell arry; same as alphas above, but with
%                  relevant xDims removed.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Oct 2022 -- Initial full revision.
%     02 Sep 2023 -- Added ability to handle alternative GP kernels.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

numIters = size(alphas,2);
for i = 1:numIters
    alphas{i} = alphas{i}(:,xDims);
    gp_params.Ds{i} = gp_params.Ds{i}(:,xDims);
    switch covType
        case {'rbf', 'exp'}
            gp_params.gams{i} = gp_params.gams{i}(xDims);
        case {'sg', 'expcos'}
            gp_params.gams{i} = gp_params.gams{i}(xDims);
            gp_params.nus{i} = gp_params.nus{i}(xDims);
    end
end
