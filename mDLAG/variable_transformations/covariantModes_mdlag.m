function [P, V, H] = covariantModes_mdlag(params, groups, varargin)
%
% [P, V, H] = covariantModes_mdlag(params, groups, ...)
%
% Description: Compute the covariant modes between a pair of groups, which
%              capture the maximal (0-lag) covariance across groups.
%                  P = V{1}'*<C{1}>*K*<C{2}>'*V{2}
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
%     P -- (xDim x xDim) array; diagonal matrix with the
%          cross-area covariance along each covariant mode.
%     V -- (1 x 2) cell array; V{i} -- (yDims(groups(i)) x xDim)
%          array; covariant modes for group groups(i), which are 
%          orthogonal but not necessarily uncorrelated.
%     H -- (1 x 2) cell array; H{i} -- (yDims(groups(i)) x xDim)
%          array; Projection of latents onto covariant modes: 
%          H{i} = V{i}'<C{i}>Ka^(0.5)
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

% Compute cross-covariance matrix
Sig_12 = params.C.means{groups(1)}*K*params.C.means{groups(2)}';

% Compute covariant modes
[V1, P, V2] = svd(Sig_12, 'econ');
P = P(1:xDim,1:xDim);
V{1} = V1(:,1:xDim);
V{2} = V2(:,1:xDim);

H{1} = V{1}'*params.C.means{groups(1)}*K^(.5);
H{2} = V{2}'*params.C.means{groups(2)}*K^(.5);
    