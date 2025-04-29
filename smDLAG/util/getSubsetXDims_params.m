function outparams = getSubsetXDims_params(inparams, xDims)
%
% outparams = getSubsetXDims_params(inparams, xDims)
%
% Description: Return a smDLAG model with only parameters corresponding to
%              the specified latents.
%
% Arguments:
%
%     inparams  -- Structure containing mDLAG model parameters. 
%                  Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         Z          -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                       array with S_j inducing point locations for
%                       latent j. Units are in time-steps.
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
%     xDims     -- (1 x numDims) array; latent state dimensions to be  
%                  retained in outparams.
%
% Outputs:
%
%     outparams -- Structure containing subset of smDLAG model parameters.
%                  Same format as inparams, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Nov 2022 -- Initial full revision.

numGroups = length(inparams.yDims);

% Initialize output structure
outparams = inparams;
outparams.xDim = length(xDims);

% State parameters
outparams.gamma = outparams.gamma(xDims);
outparams.eps = outparams.eps(xDims);
outparams.D = outparams.D(:,xDims);
outparams.Z = outparams.Z(xDims);

% ARD parameters
outparams.alpha.mean = outparams.alpha.mean(:,xDims);
outparams.alpha.b = outparams.alpha.b(:,xDims);

% Loading matrix, C
for groupIdx = 1:numGroups
    outparams.C.means{groupIdx} = outparams.C.means{groupIdx}(:,xDims);
    for yIdx = 1:outparams.yDims(groupIdx)
        if isfield(outparams.C, 'covs')
            outparams.C.covs{groupIdx}{yIdx} ...
                = outparams.C.covs{groupIdx}{yIdx}(xDims,xDims);
        end
        outparams.C.moments{groupIdx}{yIdx} ...
            = outparams.C.moments{groupIdx}{yIdx}(xDims,xDims);
    end
end
