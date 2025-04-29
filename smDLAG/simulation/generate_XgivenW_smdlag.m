function seq = generate_XgivenW_smdlag(seq, params, T, varargin)
% 
% seq = generate_XgivenW_smdlag(seq, params, T,...)
%
% Description: Generate N independent latent sequences of length T samples, 
%              given sets of inducing points.
%
% Arguments:
%
%     Required:
%
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId       -- int; unique trial (sequence) identifier
%                (inducefield) -- (xDim x 1) cell array; element j is 
%                                 a (1 x S_j) array of the inducing 
%                                 points for latent j
%
%     params  -- Structure containing mDLAG model parameters.
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
%     T         -- (1 x N) int array; T(n) gives the number of samples 
%                  (time length) for sequence n. If T is a scalar
%                  (length-1), then all sequences will be length-T.
%
%     Optional:
%
%     inducefield -- string; Name of data field in seq (default: 'wsm')
%     latentfield -- string; Name of data field in seq (default: 'xsm')
%     verbose     -- logical; Print status info (default: false)
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (inducefield) -- (xDim x 1) cell array; element j is 
%                                 a (1 x S_j) array of the inducing 
%                                 points for latent j
%                (latentfield) -- (numGroups*xDim x T) array; delayed
%                                 latent sequences
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     29 Nov 2022 -- Initial full revision.

inducefield = 'wsm';
latentfield = 'xsm';
verbose = false;
extraOpts = assignopts(who, varargin);
N = length(seq);
   
% If T is a scalar, then give all sequences the same length
if length(T) <= 1
    T = repmat(T,1,N);
end

% Group trials of same length together
Tu = unique(T);

% Initialize output structure
numGroups = length(params.yDims);
xDim = params.xDim;
for n = 1:N
    seq(n).trialId = n;
    seq(n).T = T(n);
    seq(n).(latentfield) = nan(numGroups*xDim,T(n));
end
    
% Generate all trials of the same length
for i = 1:length(Tu)
    Ti = Tu(i);
    if verbose
        fprintf('Generating all trials of length T = %d..\n', Ti);
    end
    nList = find(T == Ti);
    % Construct GP kernel matrices
    Kx = construct_Kx_smdlag(params, Ti);
    Kxw = construct_Kxw_smdlag(params, Ti);
    Kw = construct_Kw_smdlag(params);
    for n = nList
        if verbose
            fprintf('    Trial n = %d...\n', nList(n));
        end
        for j = 1:xDim
            Kxww = (Kxw{j}/Kw{j}); % (m*Ti x Sj), M -> Ti
            Kxwx = Kxww*Kxw{j}';   % (m*Ti x m*Ti), M -> Ti
            Kxwx = 0.5 * (Kxwx + Kxwx'); % Ensure symmetry
            Xj = reshape(mvnrnd((Kxww*seq(n).(inducefield){j}')', Kx{j} - Kxwx),Ti,[])'; % (numGroups x Ti)
            idx = j:xDim:xDim*numGroups;
            seq(n).(latentfield)(idx,:) = Xj;
        end
    end
end
