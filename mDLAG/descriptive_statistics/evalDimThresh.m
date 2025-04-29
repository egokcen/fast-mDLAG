function [R2, MSE] = evalDimThresh(seq, params, cutoff_sharedvar, cutoff_snr)
%
% [R2, MSE] = evalDimThresh(seq, params, cutoff_sharedvar, cutoff_snr)
%
% Description: Evaluate leave-group-out predictive performance across a
%              specified range of shared variance cutoffs. For each cutoff,
%              if a dimension does not exceed the treshold then it will
%              be zeroed out during prediction.
% 
% Arguments:
%
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- observed data
%     params  -- Structure containing mDLAG model parameters.
%     cutoff_sharedvar -- (1 x numThresh) array; list of shared variance
%                         thresholds to search over (in range [0 1]).
%     cutoff_snr       -- float; minimum signal-to-noise ratio that a
%                         group must have for ANY latents to be considered 
%                         significant
%
% Outputs:
%
%     R2      -- (1 x numThresh) array; leave-group-out R^2 values
%     MSE     -- (1 x numThresh) array; leave-group-out mean-squared errors
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Jan 2023 -- Initial full revision.

numThresh = length(cutoff_sharedvar);
numGroups = length(params.yDims);

% Initialize outputs
MSE = nan(1,numThresh);
R2 = nan(1,numThresh);
for j = 1:numThresh
    % Determine significant dimensions
    [~,sigdims,~,~] = computeDimensionalities(params, ...
                                              cutoff_sharedvar(j), ...
                                              cutoff_snr);
    % Zero-out the posterior mean and second moment of C for dimensions
    % that fall below the cutoff
    sigparams = params;
    for groupIdx = 1:numGroups
        sigparams.C.means{groupIdx}(:,~sigdims(groupIdx,:)) = 0;
        for yIdx = 1:params.yDims(groupIdx)
            sigparams.C.moments{groupIdx}{yIdx}(~sigdims(groupIdx,:),:) = 0;
            sigparams.C.moments{groupIdx}{yIdx}(:,~sigdims(groupIdx,:)) = 0;
        end
    end
    % Evaluate leave-group-out predictive performance
    [R2(j), MSE(j)] = pred_mdlag(seq, sigparams);
end