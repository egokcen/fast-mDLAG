function [R2, MSE] = pred_mdlag(seq, params)
%
% [R2, MSE] = pred_mdlag(seq, params)
%
% Description: Performs leave-group-out prediction using an existing mDLAG
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
%     MSE     -- float; leave-group-out mean-squared error
%     R2      -- float; leave-group-out R^2
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Oct 2022 -- Initial full revision.

yDims = params.yDims;
xDim = params.xDim;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

Ys_true = seq2cell2D(seq, yDims, 'datafield', 'y');

% Construct joint data matrix
Ytrue = cat(1,Ys_true{:});
[yDim, NT] = size(Ytrue);

Ys_pred = cell(1,numGroups);
for groupIdx = 1:numGroups
    Ys_pred{groupIdx} = nan(size(Ys_true{groupIdx}));
end

% Perform leave-group-out prediction
for groupIdx = 1:numGroups
    
    targetGroup = groupIdx; % Group to be left out
    sourceGroups = setdiff(1:numGroups,targetGroup); % Observed groups
    
    % Infer latent variables given source groups
    paramsSource = getSubsetGroups_params(params,sourceGroups);
    seqSource = getSubsetGroups_seq(seq,yDims,sourceGroups);
    [seqSource,~,~] = inferX(seqSource,paramsSource);
    
    % Infer latent variables for target group, given latents for source
    % groups
    seqTarget = predX(targetGroup,params,seqSource);
    Xtarget = [seqTarget.xsm];
    
    % Predict observations for target group
    paramsTarget = getSubsetGroups_params(params,targetGroup);
    Ys_pred{targetGroup} = paramsTarget.C.means{1} * Xtarget ...
        + repmat(paramsTarget.d.mean, [1 NT]);
end

% Compute performance metrics
Ypred = cat(1,Ys_pred{:});
% MSE
MSE = immse(Ypred, Ytrue);
% R2
RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
R2 = 1 - RSS / TSS;