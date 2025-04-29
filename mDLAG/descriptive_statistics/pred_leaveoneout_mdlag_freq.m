function [R2, MSE] = pred_leaveoneout_mdlag_freq(seq, params, varargin)
%
% [R2, MSE] = pred_leaveoneout_mdlag_freq(seq, params, ...)
%
% Description: Performs leave-one-out prediction using an existing mDLAG
%              model, with an approximate frequency domain approach.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- observed data
%     params  -- Structure containing mDLAG model parameters.
%
%     Optional:
%
%     eval_cutoff -- int; Evaluate performance over the time interval 
%                    (eval_cutoff+1):(T-eval_cutoff). Evaluating
%                    some middle portion of the trial can mitigate bias
%                    due to periodic edge effects. (Default: 0)     
%
% Outputs:
%
%     MSE     -- float; leave-one-out mean-squared error
%     R2      -- float; leave-one-out R^2
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     07 Jun 2024 -- Initial full revision.

eval_cutoff = 0;
assignopts(who, varargin);

% Treat each observed dimension as its own group, and then perform
% leave-group-out prediction.
params = explodeGroups_params(params);
[R2, MSE] = pred_mdlag_freq(seq, params, 'eval_cutoff', eval_cutoff);