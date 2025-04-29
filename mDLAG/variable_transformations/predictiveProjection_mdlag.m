function seq = predictiveProjection_mdlag(seq, params, groups, varargin)
%
% seq = predictiveProjection_dlag(seq, params, groups, ...)
%
% Description: Project the latent trajectories inferred by a mDLAG
%              model onto predictive modes, which capture the maximal 
%              (0-lag) predictive power between a source and target group.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     xsm (xDim x T) -- latent trajectories
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
%     groups  -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order matters: groups(1) gives the source
%                  group; groups(2) gives the target group.
%
%     Optional:
%
%     orth      -- logical; If true, project source trajectories onto an 
%                  orthogonal (as opposed to uncorrelated) basis. 
%                  (default: false)
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     seq     -- input data structure with new field
%                    xpred -- (2*xDim x T) array; trajectories projected 
%                             onto predictive modes. The first xDim latents
%                             correspond to the source group. The next xDim
%                             latents correspond to the target group.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

orth = false;
zerolag = true;
assignopts(who,varargin);

% Constants
N = length(seq);
numGroups = length(params.yDims);
xDim = params.xDim;

% Initialize output structure
for n = 1:N
    seq(n).xpred = []; 
end

% Compute predictive modes
[~, ~, U, V, ~] = predictiveModes_mdlag(params, groups, 'zerolag', zerolag);

% Project latents for each group separately
for groupIdx = 1:numGroups
    % Get appropriate latent trajectories
    groupSeq = getSubsetGroups_seq(seq, ...
                                   xDim.*ones(1,numGroups), ...
                                   groups(groupIdx), ...
                                   'datafield', 'xsm');
    X = [groupSeq.xsm];
    if orth || groupIdx > 1
        % Project target and/or source onto orthogonal basis
        Xpred = V{groupIdx}' * params.C.means{groups(groupIdx)} * X;
    else
        % Project source onto uncorrelated basis
        Xpred = U' * params.C.means{groups(groupIdx)} * X;
    end
    % Add projected latents to output structure
    seqTemp = segmentByTrial(seq,Xpred,'xpred');
    for n = 1:N
        seq(n).xpred = [seq(n).xpred; seqTemp(n).xpred];
    end
end