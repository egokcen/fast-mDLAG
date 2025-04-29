function Kxw = construct_Kxw_smdlag(params,T)
%
% Kxw = construct_Kxw_smdlag(params, T)
%
% Description: Construct GP covariance matrix between each latent and its
%              inducing points.
%
% Arguments:
%
%     params  -- Structure containing smDLAG model parameters.
%                Contains the following relevant fields:
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    Z       -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                               array with S_j inducing point locations for
%                               latent j. Units are in time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
%     T       -- int; number of timesteps
%
% Outputs:
%
%     Kxx -- (xDim x 1) cell array; Kxw{j} is a (numGroups*T) x Sj
%            GP covariance matrix between latent j and its inducing points          
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     29 Nov 2022 -- Initial full revision.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
mT           = numGroups*T;
Kxw = cell(xDim,1);
for j = 1:xDim
    Sj = length(params.Z{j});
    Kxw{j} = zeros(mT,Sj); % (M -> T) x Sj
end

for j = 1:xDim
    Tall = repmat(1:T,numGroups,1)'; % (T x m)
    Dall = repmat(params.D(:,j),1,T)'; % (T x m)
    Tdif = Tall(:) - Dall(:); % (mT x mT), M -> T
    deltaT = Tdif - params.Z{j}; % (M -> T) x Sj
    deltaTsq = deltaT.^2;
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(j)*deltaTsq);          
    end
    Kxw{j} = (1-params.eps(j))*temp; % (M -> T) x Sj
end
