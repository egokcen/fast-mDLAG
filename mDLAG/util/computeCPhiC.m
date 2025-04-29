function CPhiC = computeCPhiC(params)
%
% CPhiC = computeCPhiC(params)
%
% Description: Helper function to compute the expectation of C'*Phi*C, 
%              for each group.
%
% Arguments:
%
%     params -- Structure containing mDLAG model parameters.
%               Contains the relevant fields
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
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%
% Outputs
%
%     CPhiC  -- (1 x numGroups) cell array; CPhiC{m} is a (xDim x xDim)
%               array containing the expectation of C_m'*Phi_m*C_m
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Aug 2023 -- Initial full revision.

% Initialize relevant variables
yDims = params.yDims;
xDim = params.xDim;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

% Compute the expectation of C'*Phi*C, an intermediate term
CPhiC = cell(1,numGroups);
for groupIdx = 1:numGroups
    currGroup = block_idxs{groupIdx};
    phi_m = params.phi.mean(currGroup(1):currGroup(2));
    CPhiC{groupIdx} = zeros(xDim,xDim);
    for yIdx = 1:yDims(groupIdx)
        CPhiC{groupIdx} = CPhiC{groupIdx} + phi_m(yIdx).*params.C.moments{groupIdx}{yIdx};
    end
end
