function [R2, MSE] = pred_mdlag_freq(seq, params, varargin)
%
% [R2, MSE] = pred_mdlag_freq(seq, params, ...)
%
% Description: Performs leave-group-out prediction using an existing mDLAG
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
%     MSE     -- float; leave-group-out mean-squared error
%     R2      -- float; leave-group-out R^2
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.
%     07 Jun 2024 -- Added eval_period option.

eval_cutoff = 0;
assignopts(who, varargin);

yDims = params.yDims;
numGroups = length(yDims);

% Compute the FFT of the observed data
seq = fftseq(seq, 'y', 'yfft');

% Extract the evaluation period data
for n = 1:length(seq)
    T = seq(n).T;
    seq(n).y_eval = seq(n).y(:,(eval_cutoff+1):(T-eval_cutoff));
end
Ys_true = seq2cell2D(seq, yDims, 'datafield', 'y_eval');

% Construct joint data matrix
Ytrue = cat(1,Ys_true{:});
[~, NT] = size(Ytrue);

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
    seqSource = getSubsetGroups_seq(seq,yDims,sourceGroups,'datafield','yfft');
    [seqSource,~,~] = inferX_freq(seqSource,paramsSource);
    
    % Time-delay latents for the target group
    paramsTarget = getSubsetGroups_params(params,targetGroup);
    seqTarget = freq2time_mdlag(seqSource,paramsTarget);
    % Extract the evaluation period data
    for n = 1:length(seqTarget)
        T = seqTarget(n).T;
        seqTarget(n).xsm = seqTarget(n).xsm(:,(eval_cutoff+1):(T-eval_cutoff));
    end
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