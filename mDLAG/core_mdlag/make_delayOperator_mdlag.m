function Q = make_delayOperator_mdlag(params, f)
%
% Q = make_delayOperator_mdlag(params, f)
%
% Construct a compressed representation of the (sparse) time-delay operator
% matrix for a given frequency.
%
% Arguments:
%
%     params -- mDLAG model parameters
%     f      -- float; frequency, in units of (1/timeSteps)
%
% Outputs:
%
%     Q      -- GP time-delay operator matrix with dimensions 
%               (numGroups x xDim).
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     12 Aug 2024 -- Initial full revision.

numGroups = length(params.yDims);
Q = zeros(numGroups, params.xDim); % Time-delay operator matrix

% Fill in Q
for groupIdx = 1:numGroups
    Q(groupIdx,:) = exp(-1i*2*pi*f.*params.D(groupIdx,:));
end
