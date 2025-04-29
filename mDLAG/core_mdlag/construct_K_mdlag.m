function [K_big] = construct_K_mdlag(params,T)
%
% [K_big] = construct_K_mdlag(params, T)
%
% Description: Constructs full GP covariance matrix across all latent
%              state dimensions and timesteps.
%
% Arguments:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the following relevant fields:
% 
%                    covType -- string; type of GP covariance. The 
%                               following GP kernels are currently 
%                               supported:
%                                   'rbf'    -- Radial basis function, or 
%                                               squared exponential
%                                   'sg'     -- Spectral Gaussian
%                                   'exp'    -- Exponential
%                                   'expcos' -- Exponential-cosine
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%                    The following fields depend on GP covariance type:
%                        For 'rbf':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma)                                                    
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'sg':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma) 
%                            nu    -- (1 x xDim) array; center frequencies; 
%                                     convert to 1/time via nu./binWidth 
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'exp':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ gamma                                                    
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'expcos':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ gamma 
%                            nu    -- (1 x xDim) array; center frequencies; 
%                                     convert to 1/time via nu./binWidth 
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%
%     T       -- int; number of timesteps
%
% Outputs:
%
%     K_big   -- (xDim * numGroups * T) x (xDim * numGroups * T) array;
%                GP covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     27 Sep 2022 -- Initial full revision.
%     01 Sep 2023 -- Added spectral Gaussian compatibility.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
K_big        = zeros(xDim*numGroups*T);
mT           = numGroups*T;

Tdif = repmat(1:T,numGroups,1); % (m x T)
Tdif = repmat(Tdif(:)',mT,1) - repmat(Tdif(:),1,mT); % (mT x mT)
for j = 1:xDim
    Delayall = params.D(:,j); % (m x 1)
    Delaydif = repmat(Delayall,T,1);    % (mT x 1)
    Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT); % (mT x mT)
    deltaT = Tdif - Delaydif; 
    deltaTsq = deltaT.^2;
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(j)*deltaTsq);   
        case 'sg'
            temp = exp(-0.5*params.gamma(j)*deltaTsq).*cos(2*pi*params.nu(j)*deltaT);
        case 'exp'
            temp = exp(-params.gamma(j)*abs(deltaT));
        case 'expcos'
            temp = exp(-params.gamma(j)*abs(deltaT)).*cos(2*pi*params.nu(j)*deltaT);
    end
    K_j = (1-params.eps(j))*temp + params.eps(j)*eye(mT); % (mT x mT)
    
    idx = j:xDim:xDim*numGroups*T;
    K_big(idx,idx) = K_j;
end
