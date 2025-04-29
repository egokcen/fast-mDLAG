function Kw = construct_Kw_smdlag(params)
%
% Kw = construct_Kw_smdlag(params)
%
% Description: Constructs GP covariance matrices for each latent's set
%              of inducing points.
%
% Arguments:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the following relevant fields:
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    Z       -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                               array with S_j inducing point locations for
%                               latent j. Units are in time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
% Outputs:
%
%     Kw -- (xDim x 1) cell array; Kw{j} is a (Sj x Sj) 
%           GP covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     28 Nov 2022 -- Initial full revision.

xDim = params.xDim;
Kw = cell(xDim,1);
for j = 1:xDim
    Sj = length(params.Z{j});
    Kw{j} = zeros(Sj);
end

for j = 1:xDim
    Sj = length(params.Z{j});
    Tdif = repmat(params.Z{j}', 1, Sj) - repmat(params.Z{j}, Sj, 1);
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(j)*Tdif.^2);          
    end
    Kw{j} = (1-params.eps(j))*temp + params.eps(j)*eye(Sj); % (Sj x Sj)
end
