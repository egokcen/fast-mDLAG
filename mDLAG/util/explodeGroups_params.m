function params = explodeGroups_params(params)
%
% params = explodeGroups_params(params)
%
% Description: 'Explode' the group structure of mDLAG model parameters so 
%              that each observed dimension is its own 'group'.
%
% Arguments:
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
% Outputs:
%
%     params -- Same structure as 'params' above, but with each observed
%               dimension treated as its own group. Therefore, 
%               params.yDims = ones(1,yDim) and numGroups = yDim.               
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     06 Jun 2024 -- Initial full revision.

% Dimensionality of each group
% Keep track of old groupings
yDims_old = params.yDims;           
numGroups = length(yDims_old);

% Update groupings in params
yDim = sum(params.yDims);
yDims = ones(1,yDim);
params.yDims = yDims;

% Observed dimensions belonging to the same group will all have the same
% delay
D = [];
for groupIdx = 1:numGroups
    D = [D; repmat(params.D(groupIdx,:), [yDims_old(groupIdx), 1])];
end
params.D = D;

% Split up all rows of the loading matrix C
C_means = cell(1,yDim);
if isfield(params.C, 'covs')
    C_covs = cell(1,yDim);
end
C_moments = cell(1,yDim);
yIdx_new = 1;  % Keep track of the overall observed dimension index
for groupIdx = 1:numGroups
    for yIdx_old = 1:yDims_old(groupIdx)
        C_means{yIdx_new} = params.C.means{groupIdx}(yIdx_old,:);
        if isfield(params.C, 'covs')
            C_covs{yIdx_new} = params.C.covs{groupIdx}(yIdx_old);
        end
        C_moments{yIdx_new} = params.C.moments{groupIdx}(yIdx_old);
        yIdx_new = yIdx_new + 1;
    end
end
params.C.means = C_means;
if isfield(params.C, 'covs')
    params.C.covs = C_covs;
end
params.C.moments = C_moments;

% Observed dimensions belonging to the same group will all have the same
% ARD parameters, alpha
alpha_a = [];
alpha_b = [];
alpha_mean = [];
for groupIdx = 1:numGroups
    alpha_a = [alpha_a; repmat(params.alpha.a(groupIdx), [yDims_old(groupIdx), 1])];
    alpha_b = [alpha_b; repmat(params.alpha.b(groupIdx,:), [yDims_old(groupIdx), 1])];
    alpha_mean = [alpha_mean; repmat(params.alpha.mean(groupIdx,:), [yDims_old(groupIdx), 1])];
end
params.alpha.a = alpha_a;
params.alpha.b = alpha_b;
params.alpha.mean = alpha_mean;
