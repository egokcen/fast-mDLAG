function [P, targetVar, U, V, H] = predictiveModes_mdlag(params, groups, varargin)
%
% [P, targetVar, U, V, H] = predictiveModes_mdlag(params, groups, ...)
%
% Description: Compute the predictive modes between a source and target
%              group, which capture maximal (0-lag) predictive power from
%              source to target.
%                  P = V{1}'*Sig_11^(-0.5)*<C{1}>*K*<C{2}>'*V{2}
%                    = U{1}'*<C{1}>*K*<C{2}>'*V{2}
%                    = H{1}*H{2}'
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
%     groups  -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order matters: groups(1) gives the source
%                  group; groups(2) gives the target group.
%
%     Optional:
%
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     P  -- (xDim x xDim) array; diagonal matrix with the
%           cross-area predictive power along each predictive mode.
%     targetVar -- float; total variance of the target group.
%     U  -- (sourceDim x xDim) array;
%           predictive modes for the source group, which are uncorrelated 
%           but not necessarily orthogonal. U = Sig_11^(-.5)*V{1}
%     V  -- (1 x 2) cell array; V{i} -- (yDims(groups(i)) x xDim)
%           array; predictive modes for group groups(i), which are 
%           orthogonal but not necessarily uncorrelated.
%     H  -- (1 x 2) cell array; H{i} -- (yDims(groups(i)) x xDim)
%           array; Projection of latents onto predictive modes: 
%           H{1} = U'<C{1}>K; H{2} = V{2}'<C{2}>K
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

zerolag = true;
assignopts(who,varargin);
numGroups = length(groups);

% Initialize outputs
P = [];
targetVar = nan;
U = [];
V = cell(1,numGroups);
H = cell(1,numGroups);

xDim = params.xDim;
if zerolag
    % Get the 0-lag portion of the GP kernel matrix, K(t,t)
    K = getZeroLagK(params,groups);
else
    % Otherwise, set K to the identity matrix, thereby maintaining
    % time delays in the computation of modes.
    K = eye(xDim);
end

% Compute least-squares matrix
params1 = getSubsetGroups_params(params,groups(1));
Sig_11 = params1.C.means{1}*params1.C.means{1}' +  diag(params1.phi.mean);
params2 = getSubsetGroups_params(params,groups(2));
Sig_22 = params2.C.means{1}*params2.C.means{1}' +  diag(params2.phi.mean);
Sig_12 = params1.C.means{1}*K*params2.C.means{1}';
Corr_12 = Sig_11^(-0.5) * Sig_12;

% Compute predictive modes
[V1, P, V2] = svd(Corr_12, 'econ');
P = P(1:xDim,1:xDim);
targetVar = trace(Sig_22);
V{1} = V1(:,1:xDim);
V{2} = V2(:,1:xDim);

U = Sig_11^(-0.5)*V{1};

H{1} = U'*params.C.means{groups(1)}*K^(.5);
H{2} = V{2}'*params.C.means{groups(2)}*K^(.5);
    