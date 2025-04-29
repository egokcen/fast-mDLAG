function precomp = makePrecomp_GP(seq,params)
%
% precomp = makePrecomp_GP(seq,params)
%
% Description: Precompute several quantities that get re-used each GP
%              gradient update.
%
% Arguments:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%                     wsm          -- (xDim x 1) cell array; element j is a
%                                     (1 x S_j) array of the inducing 
%                                     points for latent j
%
%     params -- Structure containing mDLAG model parameters.
%               Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
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
%         W          -- data structure whose jth entry, corresponding
%                       to a group of trials of the same length, has fields
%                            T      -- int; number of time steps for this
%                                      trial group
%                            cov    -- (S x S) posterior covariance of the  
%                                      inducing points for all latents.
%                            moment -- (S x S) posterior second moment of 
%                                      the inducing points for all latents.
%
% Outputs
%     precomp -- Structure whose ith entry contains precomputations for 
%                the i-th latent state:
%                xIdx      -- int; index of i-th latent state
%                numGroups -- int; number of observed groups
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        nList -- List of trial IDs belonging to this group
%                        T     -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        KxwKw_inv -- (1 x xDim) cell array; Product
%                                     of each Kxw and inv(Kw)
%                        Wmoment -- (1 x xDim) cell array; Precomputed 
%                                   posterior second moments of W that
%                                   involve latent i
%                        CPhiYW -- (numGroups*T x Si) array; intermediate
%                                  term
%                        CPhiC -- (1 x xDim) cell array; each element
%                                 is a (numGroups x 1) array; the 
%                                 diagonal elements of a restructured CPhiC
%                        CPhiC_big -- (1 x xDim) cell array; each element
%                                     is a (numGroups*T x 1) array;
%                                     contains copies of the diagonal 
%                                     elements of CPhiC
%                        WWKwKwx -- (1 x xDim) cell array; each
%                                   element is a (Si x numGroups*T)
%                                   array; another intermediate term
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Dec 2022 -- Initial full revision.
%     11 Aug 2024 -- Overhauled with learnGPparams, grad_GPparams to revise
%                    several supoptimal computations.

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
yDim = sum(yDims);
numGroups = length(yDims);
obs_block_idxs = get_block_idxs(yDims);
S = nan(1,xDim);
for i = 1:xDim
    S(i) = length(params.Z{i}); 
end
Z_block_idxs = get_block_idxs(S); % Index blocks of inducing points

% Assign some helpful precomp items
for i = 1:xDim
    precomp(i).xIdx = i;
    precomp(i).numGroups = numGroups;
end
% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);
% Loop once for each state dimension (each GP)
for i = 1:xDim
    for j = 1:length(Tu)
        T = Tu(j);
        precomp(i).Tu(j).nList = find(Tall == T);
        precomp(i).Tu(j).T = T;
        precomp(i).Tu(j).numTrials = length(precomp(i).Tu(j).nList);
        precomp(i).Tu(j).Wmoment  = cell(1,xDim);
        precomp(i).Tu(j).CPhiYW  = zeros(numGroups*T,S(i));
        precomp(i).Tu(j).CPhiC = cell(1,xDim);
        precomp(i).Tu(j).CPhiC_big = cell(1,xDim);
        precomp(i).Tu(j).WWKwKwx = cell(1,xDim);
    end
end

% Precompute products between GP covariances
for j = 1:length(Tu)
    T = Tu(j);
    Kxw = construct_Kxw_smdlag(params, T);
    Kw = construct_Kw_smdlag(params);
    precompK.Tu(j).KxwKw_inv = cell(1,xDim);
    for xIdx = 1:xDim
        try
            precompK.Tu(j).KxwKw_inv{xIdx} ...
                = Kxw{xIdx} * invChol_mex(Kw{xIdx});
        catch
            precompK.Tu(j).KxwKw_inv{xIdx} ...
                = Kxw{xIdx} / Kw{xIdx};
        end
    end
    
    % Add these to the primary precomp structure
    for i = 1:xDim
        precomp(i).Tu(j).KxwKw_inv = precompK.Tu(j).KxwKw_inv;
    end
end

% Fill out Wmoment
% Loop once for each state dimension (each GP)
for i1 = 1:xDim
    currLatent1 = Z_block_idxs{i1};
    % Loop once for each trial length (each of Tu)
    for j = 1:length(Tu)
        for i2 = 1:xDim
            currLatent2 = Z_block_idxs{i2};
            precomp(i1).Tu(j).Wmoment{i2} ...
                = params.W(j).moment(currLatent1(1):currLatent1(2), ...
                                     currLatent2(1):currLatent2(2));
        end
    end
end

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
% Restructure CPhiC so that it matches the organization of W
tmp = cell(xDim,xDim);
for k = 1:xDim
    for j = 1:xDim
        % Submatrices are diagonal, so keep only the diagonal elements
        tmp{k,j} = zeros(numGroups,1);
        for groupIdx = 1:numGroups
            tmp{k,j}(groupIdx) = CPhiC{groupIdx}(k,j); 
        end
    end
end
CPhiC = tmp;

for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T); 
    
    % An intermediate term
    CPhiY0 = zeros(numGroups*T,length(nList),xDim); % (mT x N x p), M -> T
    group_block_idxs = get_block_idxs(repmat(T,numGroups));
    for groupIdx = 1:numGroups
        currObsGroup = obs_block_idxs{groupIdx}; % Observation groups
        groupBlock = group_block_idxs{groupIdx}(1):group_block_idxs{groupIdx}(2);
        for xIdx = 1:xDim
            for i = 1:length(nList)
                n = nList(i);
                CPhiY0(groupBlock,i,xIdx) = ...
                    (seq(n).y(currObsGroup(1):currObsGroup(2),:) ...
                    - params.d.mean(currObsGroup(1):currObsGroup(2)))' ...
                    * CPhi{groupIdx}(xIdx,:)';
            end
        end
    end
    
    % CPhiYW and WWKwKwx
    for xIdx1 = 1:xDim
        wsmMat = zeros(S(xIdx1),length(nList));
        for i = 1:length(nList)
            n = nList(i);
            wsmMat(:,i) = seq(n).wsm{xIdx1}'; % (Sj x length(nList))
        end
        
        precomp(xIdx1).Tu(j).CPhiYW = CPhiY0(:,:,xIdx1) * wsmMat';
        
        for xIdx2 = 1:xDim
            precomp(xIdx1).Tu(j).CPhiC{xIdx2} = CPhiC{xIdx1,xIdx2}; % (M x 1)
            precomp(xIdx1).Tu(j).CPhiC_big{xIdx2} = kron(CPhiC{xIdx1,xIdx2},ones(T,1)); % (MT x 1), M -> T
            precomp(xIdx1).Tu(j).WWKwKwx{xIdx2} ...
                = precomp(xIdx1).Tu(j).Wmoment{xIdx2} ...
                * precomp(xIdx1).Tu(j).KxwKw_inv{xIdx2}';
        end
    end
end
