function [seq, params, lbterm, Xmoment] = inferXW(seq, params, varargin)
%
% [seq, params, lbterm, Xmoment] = inferXW(seq, params, ...)
%
% Description: Infer inducing points W and latents X given observations Y 
%              (in seq) and smDLAG model parameters in params.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- neural data
%
%     params  -- Structure containing smDLAG model parameters. 
%                Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         Z          -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                       array with S_j inducing point locations for
%                       latent j. Units are in time-steps.
%         d.mean     -- (yDim x 1) array; posterior mean of mean parameter
%         d.cov      -- (yDim x 1) array; diagonal elements of the
%                       posterior covariance matrix of d
%         C.means    -- (numGroups x 1) cell array; yDims(groupIdx) x xDim
%                       mean loadings matrix for each group
%         C.covs     -- (numGroups x 1) cell array; C.covs{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       covariance of a row of C.
%         C.moments  -- (numGroups x 1) cell array; C.moments{groupIdx) is a
%                       (yDims(groupIdx) x 1) cell array, and each element
%                       is a (xDim x xDim) matrix giving the posterior
%                       second moment of a row of C.
%         alpha.a    -- (numGroups x 1) array; shape parameters of alpha 
%                       posterior
%         alpha.b    -- (numGroups x xDim) array; scale parameters of 
%                       alpha posterior
%         alpha.mean -- (numGroups x xDim) array; mean precisions of
%                       loading weights (for ARD); alpha.a ./ alpha.b
%         phi.a      -- float; shape parameter of phi posterior
%         phi.b      -- (yDim x 1) array; scale parameters of phi posterior
%         phi.mean   -- (yDim x 1) array; mean precisions of observations; 
%                       alpha.a ./ alpha.b
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%     
%     Optional:
%
%     inferX -- logical; Set true to infer X in addition to W, false
%               to only infer W (default: true)
%
% Outputs:
%
%     seq    -- data structure whose nth entry has new fields (these fields 
%               are added to existing fields in the seq input argument)
%                 (if inferX is true)
%                 xsm  -- ((numGroups*xDim) x T) array; posterior mean 
%                         at each timepoint
%                 wsm  -- (xDim x 1) cell array; element j is a (1 x S_j)
%                         array of the inducing points for latent j
%     params -- Structure containing smDLAG model parameters with new
%               fields
%                 NOTE: For smDLAG, posterior covariances of X are the same 
%                       for trials of the same length.
%                 W     -- data structure whose jth entry, corresponding
%                          to a group of trials of the same length, has
%                          fields
%                            T      -- int; number of time steps for this
%                                      trial group
%                            cov    -- (S x S) posterior covariance of the  
%                                      inducing points for all latents.
%                            moment -- (S x S) posterior second moment of 
%                                      the inducing points for all latents.
%     lbterm -- float; A term that contributes to the variational lower 
%               bound
%     (if inferX is true)
%     Xmoment -- (1 x numGroups) cell array; Xmoment{m} is a (xDim x xDim)
%                array containing second moments of X for the mth group. 
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Nov 2022 -- Initial full revision.
%     10 Aug 2024 -- Overhauled several suboptimal computations.

inferX = true;
assignopts(who, varargin);

% Initialize other relevant variables
yDims = params.yDims;  
yDim = sum(yDims);
xDim = params.xDim;
numGroups = length(yDims);
obs_block_idxs = get_block_idxs(yDims); % Index blocks of observed variables
S = nan(1,xDim);
for j = 1:xDim
    S(j) = length(params.Z{j}); 
end
Z_block_idxs = get_block_idxs(S); % Index blocks of inducing points
Xmoment = cell(1,numGroups);
if inferX
    for groupIdx = 1:numGroups
        Xmoment{groupIdx} = zeros(xDim);
    end
end
Kxw_mt = zeros(xDim,sum(S));  % Preallocate this array for later. % (p x S), p -> S

% Compute C'*Phi and the expectation of C'*Phi*C, intermediate terms
CPhi = cell(1,numGroups);
CPhiC = cell(1,numGroups);
for groupIdx = 1:numGroups
    currGroup = obs_block_idxs{groupIdx};
    phi_m = params.phi.mean(currGroup(1):currGroup(2));
    CPhi{groupIdx} = params.C.means{groupIdx}' .* phi_m';
    CPhiC{groupIdx} = zeros(xDim,xDim);
    for yIdx = 1:yDims(groupIdx)
        CPhiC{groupIdx} = CPhiC{groupIdx} + phi_m(yIdx).*params.C.moments{groupIdx}{yIdx};
    end
end

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.
lbterm = 0;
for j = 1:length(Tu)
    T = Tu(j);
    params.W(j).T = T;
    % Process all trials with length T
    nList = find(Tall == T);

    % Subtract estimated mean from observations
    Y0 = zeros(yDim,length(nList),T); % (q x N x T)
    for i = 1:length(nList)
        n = nList(i);
        Y0(:,i,:) = permute(seq(n).y - params.d.mean, [1 3 2]);
    end
    
    % Construct GP kernel matrices
    Kxw = construct_Kxw_smdlag(params, T);
    Kw = construct_Kw_smdlag(params);
    Kw_big = blkdiag(Kw{:});

    group_block_idxs = get_block_idxs(repmat(T,numGroups)); % Index blocks of Kxw
    
    % Initialize the posterior mean and covariance of W
    SigW = Kw_big;
    wsmMat = zeros(sum(S), length(nList));
    for groupIdx = 1:numGroups
        currObsGroup = obs_block_idxs{groupIdx}; % Observation groups
        groupBlock = group_block_idxs{groupIdx}(1):group_block_idxs{groupIdx}(2);
        for t = 1:T
            % Construct submatrices of Kxw corresponding to the current group
            % and time point
            for xIdx = 1:xDim
                Zblock = Z_block_idxs{xIdx};
                Kxw_mt(xIdx,Zblock(1):Zblock(2)) = Kxw{xIdx}(groupBlock(t),:); % t-th row of Kxw_m
            end
            % Compute an intermediate term of the posterior covariance of W
            SigW = SigW + Kxw_mt' * CPhiC{groupIdx} * Kxw_mt;
            % Compute an intermediate term of the posterior mean of W
            wsmMat = wsmMat ...
                + Kxw_mt' * CPhi{groupIdx} * Y0(currObsGroup(1):currObsGroup(2),:,t);
        end
    end
    % Intermediate term -- saves computation between mean and covariance
    try
        SigW_Kinv = Kw_big * invChol_mex(SigW);
    catch
        SigW_Kinv = Kw_big / SigW;
    end
    % Finish up posterior covariance of W
    SigW = SigW_Kinv * Kw_big; % (S x S)
    SigW = 0.5 * (SigW + SigW');
    params.W(j).cov = SigW;

    % Finish up posterior mean of W
    wsmMat = SigW_Kinv * wsmMat; % S x length(nList)
    % Restructure W and collect it in seq
    for i = 1:length(nList)
        n = nList(i);
        seq(n).wsm = cell(xDim,1);
        for xIdx = 1:xDim
            currLatent = Z_block_idxs{xIdx};
            seq(n).wsm{xIdx} = wsmMat(currLatent(1):currLatent(2),i)';
        end
    end
    % Posterior second moment of W
    params.W(j).moment = length(nList).*SigW + wsmMat * wsmMat'; % (S x S)
    
    % Lower bound term
    lbterm = lbterm + length(nList)*logdet(SigW);
    logdetK = 0; % Intermediate terms
    trKWW = 0;
    Kinv = cell(1,xDim);
    for xIdx = 1:xDim
        currLatent = Z_block_idxs{xIdx};
        logdetK = logdetK + length(nList)*logdet(Kw{xIdx});
        try
            Kinv{xIdx} = invChol_mex(Kw{xIdx});
        catch
            Kinv{xIdx} = inv(Kw{xIdx});
        end
        trKWW = trKWW + Kinv{xIdx}(:)' ...
            * reshape(params.W(j).moment(currLatent(1):currLatent(2),currLatent(1):currLatent(2)),[],1);
    end
    Kinv_big = blkdiag(Kinv{:});
    lbterm = lbterm - (logdetK + trKWW);
    
    if inferX
        % We need not compute the entire posterior covariance/moment of X,
        % only covariance between latents and groups at each time point
        for groupIdx = 1:numGroups
            groupBlock = group_block_idxs{groupIdx}(1): group_block_idxs{groupIdx}(2);
            for t = 1:T
                % Construct submatrices of Kxw corresponding to the current group
                % and time point
                for xIdx = 1:xDim
                    Zblock = Z_block_idxs{xIdx};
                    Kxw_mt(xIdx,Zblock(1):Zblock(2)) = Kxw{xIdx}(groupBlock(t),:); % t-th row of Kxw_m
                end
                Kxw_Kinv = Kxw_mt * Kinv_big;  % (p x S)
                Xmoment{groupIdx} = Xmoment{groupIdx} ...
                    + length(nList).*(eye(xDim) - Kxw_Kinv * Kxw_mt') ...
                    + Kxw_Kinv * params.W(j).moment * Kxw_Kinv';
            end
        end

        % Compute posterior mean of X
        % Restructure X to match the structure expected in other functions
        for i = 1:length(nList)
            n = nList(i);
            seq(n).xsm = zeros(numGroups*xDim,T);
            for xIdx = 1:xDim
                currLatent = Z_block_idxs{xIdx};
                xsm_j = Kxw{xIdx} * (Kinv{xIdx} * wsmMat(currLatent(1):currLatent(2),i)); % (MT x 1), M -> T
                for groupIdx = 1:numGroups
                    currGroup = group_block_idxs{groupIdx};
                    xsm_jm = xsm_j(currGroup(1):currGroup(2)); % (T x 1)
                    seq(n).xsm((groupIdx-1)*xDim+xIdx,:) = xsm_jm'; % (Mp x 1), M -> p
                end
            end
        end
    end
end
