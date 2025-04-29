function [Ds,gams,Zs,alphas] = getSubsetXDims_trackedParams(Ds,gams,Zs,alphas,xDims)
%
% [Ds,gams,Zs,alphas] = getSubsetXDims_trackedParams(Ds,gams,Zs,alphas,xDims)
%
% Description: Remove unwanted latents from smDLAG model parameters that
%              are tracked throughout fitting.
%
% Arguments:
%
%     Ds        -- (1 x numIters) cell array; the estimated delay matrix
%                  (D) after each EM iteration.
%     gams      -- (1 x numIters) cell arry; estimated gamma after each EM
%                  iteration.
%     Zs        -- (1 x numIters) cell array; estimated inducing point
%                  locations after each EM iteration.
%     alphas    -- (1 x numIters) cell arry; estimated ARD parameters
%                  (alpha.mean) after each EM iteration.
%     xDims     -- (1 x numDims) array; latent state dimensions to be  
%                  retained in outparams.
%
% Outputs:
%
%     Ds        -- (1 x numIters) cell array; same as Ds above, but with
%                  relevant xDims removed.
%     gams      -- (1 x numIters) cell arry; same as gams above, but with
%                  relevant xDims removed.
%     Zs        -- (1 x numIters) cell array; same as Zs above, but with
%                  relevant xDims removed.
%     alphas    -- (1 x numIters) cell arry; same as alphas above, but with
%                  relevant xDims removed.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Dec 2022 -- Initial full revision.

numIters = length(Ds);
for i = 1:numIters
    Ds{i} = Ds{i}(:,xDims);
    gams{i} = gams{i}(xDims);
    Zs{i} = Zs{i}(xDims);
    alphas{i} = alphas{i}(:,xDims);
end
