function Q = make_delayOperator_full_mdlag(params, f)
%
% Q = make_delayOperator_full_mdlag(params, f)
%
% Construct the full time-delay operator matrix for a given frequency.
%
% Arguments:
%
%     params -- mDLAG model parameters
%     f      -- float; frequency, in units of (1/timeSteps)
%
% Outputs:
%
%     Q      -- GP time-delay operator matrix with dimensions 
%               (numGroups*xDim x xDim).
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     23 Aug 2023 -- Initial full revision.

numGroups = length(params.yDims);
Q = cell(1,numGroups); % Time-delay operator matrix for each group

% Fill in Q
for groupIdx = 1:numGroups
    Q{groupIdx} = diag(exp(-1i*2*pi*f.*params.D(groupIdx,:)).');
end

% Collect time-delay operator matrices across groups
Q = vertcat(Q{:});
