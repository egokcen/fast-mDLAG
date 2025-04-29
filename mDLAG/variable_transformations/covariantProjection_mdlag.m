function seq = covariantProjection_mdlag(seq, params, groups, varargin)
%
% seq = covariantProjection_mdlag(seq, params, groups, ...)
%
% Description: Project the latent trajectories inferred by a mDLAG
%              model onto covariant modes, which capture the maximal 
%              (0-lag) covariance between a pair of groups.
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
%                analyze. Order does not change the computation, but it
%                does change the order of groups in the outputs.
%
%     Optional:
%
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     seq     -- input data structure with new field
%                xcov  -- (2*xDim x T) array; trajectories projected onto 
%                         covariant modes. The first xDim latents 
%                         correspond to group groups(1). The next xDim 
%                         latents correspond to group groups(2).
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

zerolag = true;
assignopts(who,varargin);

% Constants
N = length(seq);
numGroups = length(params.yDims);
xDim = params.xDim;

% Initialize output structure
for n = 1:N
    seq(n).xcov = []; 
end

% Compute covariant modes
[~, V, ~] = covariantModes_mdlag(params, groups, 'zerolag', zerolag);

% Project latents for each group separately
for groupIdx = 1:length(groups)
    % Get appropriate latent trajectories
    groupSeq = getSubsetGroups_seq(seq, ...
                                   xDim.*ones(1,numGroups), ...
                                   groups(groupIdx), ...
                                   'datafield', 'xsm');
    X = [groupSeq.xsm];
    % Project onto orthogonal basis
    Xcov = V{groupIdx}' * params.C.means{groups(groupIdx)} * X;
    % Add projected latents to output structure
    seqTemp = segmentByTrial(seq,Xcov,'xcov');
    for n = 1:N
        seq(n).xcov = [seq(n).xcov; seqTemp(n).xcov];
    end
end
