function [R2, MSE] = pred_leaveoneout_mdlag(seq, params)
%
% [R2, MSE] = pred_leaveoneout_mdlag(seq, params)
%
% Description: Performs leave-one-out prediction using an existing mDLAG
%              model.
%
% Arguments:
%
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- observed data
%     params  -- Structure containing mDLAG model parameters.
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
%     06 Jun 2024 -- Initial full revision.

% Treat each observed dimension as its own group, and then perform
% leave-group-out prediction.
params = explodeGroups_params(params);
[R2, MSE] = pred_mdlag(seq, params);
