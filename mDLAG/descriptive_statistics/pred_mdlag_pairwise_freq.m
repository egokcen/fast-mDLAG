function [R2, MSE] = pred_mdlag_pairwise_freq(seq, params)
%
% [R2, MSE] = pred_mdlag_pairwise_freq(seq, params)
%
% Description: Performs prediction between each pair of groups using an 
%              existing mDLAG model, with an approximate frequency domain
%              approach.
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
%     R2.uni  -- (numGroups x numGroups) array; R2(i,j) gives the R^2 
%                 value when predicting group j from group i
%                 (unidirectional prediction, asymmetric)
%     R2.agg  -- (numGroups x numGroups) array; R2(i,j) gives the aggregate
%                R^2 value when predicting between groups i and j
%                (symmetric; only upper triangular portion is filled in)
%     MSE.uni -- (numGroups x numGroups) array; MSE(i,j) gives the
%                mean-squared error when predicting group j from group i
%                (unidirectional prediction, asymmetric)
%     MSE.agg -- (numGroups x numGroups) array; MSE(i,j) gives the aggregate
%                MSE value when predicting between groups i and j
%                (symmetric; only upper triangular portion is filled in)
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

yDims = params.yDims;
numGroups = length(yDims);
pairs = nchoosek(1:numGroups,2);
numPairs = size(pairs,1);

Ys_true = seq2cell2D(seq, yDims, 'datafield', 'y');
seq = fftseq(seq, 'y', 'yfft');

% Construct joint data matrix
[~, NT] = size(Ys_true{1});

% Perform pairwise prediction
R2.uni = nan(numGroups,numGroups);
R2.agg = nan(numGroups,numGroups);
MSE.agg = nan(numGroups,numGroups);
MSE.uni = nan(numGroups,numGroups);

for pairIdx = 1:numPairs
    
    pair = pairs(pairIdx,:);
    Ys_pred = cell(1,length(pair));
    
    for targetGroup = pair
        
        sourceGroup = setdiff(pair, targetGroup); % Observed groups
        
        % Infer latent variables given source group
        paramsSource = getSubsetGroups_params(params,sourceGroup);
        seqSource = getSubsetGroups_seq(seq,yDims,sourceGroup,'datafield','yfft');
        [seqSource,~,~] = inferX_freq(seqSource,paramsSource);

        % Time-delay latents for the target group
        paramsTarget = getSubsetGroups_params(params,targetGroup);
        seqTarget = freq2time_mdlag(seqSource,paramsTarget);
        Xtarget = [seqTarget.xsm];

        % Predict observations for target group
        paramsTarget = getSubsetGroups_params(params,targetGroup);
        Ys_pred{targetGroup} = paramsTarget.C.means{1} * Xtarget ...
            + repmat(paramsTarget.d.mean, [1 NT]);

        % Compute performance metrics for this target group
        % MSE
        MSE.uni(sourceGroup,targetGroup) = immse(Ys_pred{targetGroup}, Ys_true{targetGroup});
        % R2
        RSS = sum( sum( ( Ys_true{targetGroup} - Ys_pred{targetGroup} ).^2, 1 ) );
        TSS = sum( sum( ( Ys_true{targetGroup} - repmat( mean(Ys_true{targetGroup},2), [1 size(Ys_true{targetGroup},2)] ) ).^2, 1 ) );
        R2.uni(sourceGroup,targetGroup) = 1 - RSS / TSS;
        
    end
    
    % Compute aggregate performance metrics
    Ypred = cat(1,Ys_pred{:});
    Ytrue = cat(1,Ys_true{pair});
    % MSE
    MSE.agg(pair(1),pair(2)) = immse(Ypred, Ytrue);
    % R2
    RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
    TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
    R2.agg(pair(1),pair(2)) = 1 - RSS / TSS;
    
end