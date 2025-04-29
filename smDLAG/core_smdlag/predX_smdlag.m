function seqPred = predX_smdlag(predGroup, predParams, seqSource)
%
% seqPred = predX_smdlag(predGroup, predParams, seqSource)
%
% Description: Infer latent variables for the group 'predGroup' given 
%              the latent variables for the remaining groups in seqSource.
%
% Arguments:
%
%     predGroup  -- int; index of group to be predicted
%     predParams -- Structure containing smDLAG model parameters of the
%                   target group. Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (1 x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         Z          -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                       array with S_j inducing point locations for
%                       latent j. Units are in time-steps.
%         xDim       -- int; number of latent variables
%     seqSource  -- data structure, whose nth entry (corresponding to
%                   the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     wsm -- (xDim x 1) cell array; element j is a 
%                            (1 x S_j) array of the inducing points for 
%                            latent j
%
% Outputs:
%
%     seqPred -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId        -- unique trial identifier
%                    T (1 x 1)      -- number of timesteps
%                    xsm (xDim x T) -- latent variables for predGroup
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Dec 2022 -- Initial full revision.

% Initialize relevant variables
xDim = predParams.xDim;

% Initialize output structure
for n = length(seqSource)
    seqPred(n).trialId = seqSource(n).trialId;
    seqPred(n).T = seqSource(n).T;
    seqPred(n).xsm = zeros(xDim,seqSource(n).T);
end

% Group trials of same length together
Tall = [seqSource.T];
Tu = unique(Tall);

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.

for j = 1:length(Tu)
    T = Tu(j);
    
    % Construct GP kernel matrices
    Kxw = construct_Kxw_smdlag(predParams, T);
    Kw = construct_Kw_smdlag(predParams);
    
    % Process all trials with length T
    nList    = find(Tall == T);

    % Predict latent variables for predGroup from the latent variables in 
    % sourceGroups
    for n = nList
        for xIdx = 1:xDim
            seqPred(n).xsm(xIdx,:) = (Kxw{xIdx} / Kw{xIdx}) * seqSource(n).wsm{xIdx}';
        end
    end
    
end
