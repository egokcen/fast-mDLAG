function seq = generate_W_smdlag(params, N, varargin)
% 
% seq = generate_W_smdlag(params, N,...)
%
% Description: For each latent j, generate N independent inducing point 
%              sequences of length S_j samples, according to a zero-mean 
%              Gaussian Process defined by the smDLAG model.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing smDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    Cs      -- (1 x numGroups) cell array; List of factor 
%                               loadings 
%                               {(y1Dim x xDim), (y2Dim x xDim), ...}
%                    alphas  -- (numGroups x xDim) array; alpha parameter 
%                               values
%                    phis    -- (1 x numGroups) cell array; List of 
%                               observation precision parameters 
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    ds      -- (1 x numGroups) cell array; List of data 
%                               means
%                               {(y1Dim x 1), (y2Dim x 1), ...}
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
%     N         -- int; number of sequences
%
%     Optional:
%
%     inducefield -- string; Name of data field in seq (default: 'wsm')
%     verbose     -- logical; Print status info (default: false)
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId       -- int; unique trial (sequence) identifier
%                (inducefield) -- (xDim x 1) cell array; element j is 
%                                 a (1 x S_j) array of the inducing 
%                                 points for latent j
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     28 Nov 2022 -- Initial full revision.

inducefield = 'wsm';
verbose = false;
extraOpts = assignopts(who, varargin);

% Initialize output structure
xDim = params.xDim;
for n = 1:N
    seq(n).trialId = n;
    seq(n).(inducefield) = cell(xDim,1);
    for j = 1:xDim
        Sj = length(params.Z{j});
        seq(n).(inducefield){j} = nan(1,Sj); 
    end
end

% Generate GP kernel matrix
Kw = construct_Kw_smdlag(params);
    
% Generate trials for each latent
for j = 1:xDim
    Sj = length(params.Z{j});
    if verbose
        fprintf('Generating all trials of latent j = %d..\n', j);
    end
    for n = 1:N
        if verbose
            fprintf('    Trial n = %d...\n', n);
        end
        seq(n).(inducefield){j} = mvnrnd(zeros(1,Sj), Kw{j});
    end
end
