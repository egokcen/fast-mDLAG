function [seq, params, lbterm] = inferXW_withXcov(seq, params, varargin)
%
% [seq, params, lbterm] = inferXW_withXcov(seq, params, ...)
%
% Description: Infer inducing points W and latents X given observations Y 
%              (in seq) and smDLAG model parameters in params. Compute also
%              the posterior covariance of X (not needed to fit an smDLAG
%              model).
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
%                 (if inferX is true)
%                 X.cov -- data structure whose jth entry, corresponding to
%                          a group of trials of the same length, has 
%                          fields
%                            T     -- int; number of time steps for this
%                                     trial group
%                            Vsm   -- (xDim*numGroups x xDim*numGroups x T)
%                                     array; posterior covariance at each 
%                                     timepoint
%                            VsmGP -- (numGroups*T x numGroups*T x xDim) 
%                                     array; posterior covariance of each 
%                                     GP
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
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Nov 2022 -- Initial full revision.

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

% Total number of latents, for all groups
mp = numGroups * xDim;

% Restructure C so that it matches the organization of W
tmp = cell(1,xDim);
for j = 1:xDim
    tmp{j} = cell(1,numGroups);
    for groupIdx = 1:numGroups
        tmp{j}{groupIdx} = params.C.means{groupIdx}(:,j); 
    end
    tmp{j} = blkdiag(tmp{j}{:});
end
C = cell2mat(tmp); % (yDim x numGroups*xDim), q x (p -> M)
% Compute C'*Phi, an intermediate term
CPhi = C' * diag(params.phi.mean); % (numGroups*xDim x yDim), (p -> M) x (M -> q_m)

% Compute the expectation of C'*Phi*C, an intermediate term
CPhiC = cell(1,numGroups);
for groupIdx = 1:numGroups
    currGroup = obs_block_idxs{groupIdx};
    phi_m = params.phi.mean(currGroup(1):currGroup(2));
    CPhiC{groupIdx} = zeros(xDim,xDim);
    for yIdx = 1:yDims(groupIdx)
        CPhiC{groupIdx} = CPhiC{groupIdx} + phi_m(yIdx).*params.C.moments{groupIdx}{yIdx};
    end
end
% Restructure CPhiC so that it matches the organization of W
tmp = cell(xDim,xDim);
for k = 1:xDim
    for j = 1:xDim
        tmp{k,j} = zeros(numGroups);
        for groupIdx = 1:numGroups
            tmp{k,j}(groupIdx,groupIdx) = CPhiC{groupIdx}(k,j); 
        end
    end
end
CPhiC = cell2mat(tmp); % (numGroups*xDim x numGroups*xDim), p -> M

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
    
    % Construct GP kernel matrices
    Kxw = construct_Kxw_smdlag(params, T);
    Kw = construct_Kw_smdlag(params);
    
    % Compute the full posterior covariance of W
    CPhiC_big = kron(CPhiC, eye(T)); % (MpT x MpT), p -> M -> T
    % Intermediate term -- saves computation between mean and covariance
    try
        SigW_Kinv = blkdiag(Kw{:}) * invChol_mex(blkdiag(Kw{:}) + blkdiag(Kxw{:})' * CPhiC_big * blkdiag(Kxw{:}));
    catch
        SigW_Kinv = blkdiag(Kw{:}) / (blkdiag(Kw{:}) + blkdiag(Kxw{:})' * CPhiC_big * blkdiag(Kxw{:}));
    end
    SigW = SigW_Kinv * blkdiag(Kw{:}); % (S x S)
    SigW = 0.5 * (SigW + SigW');
    params.W(j).cov = SigW;
    
    % Process all trials with length T
    nList    = find(Tall == T);
    
    % Subtract estimated mean from observations and structure to match the
    % structure of CPhi
    Y0 = zeros(yDim*T,length(nList));
    for i = 1:length(nList)
        n = nList(i);
        Y0(:,i) = reshape((seq(n).y - params.d.mean)',[],1); % (qT x 1), M -> q_m -> T
    end

    % Posterior mean of W
    CPhi_big = kron(CPhi, eye(T)); % (MpT x qT), (p -> M -> T) x (M -> q_m -> T)
    wsmMat = SigW_Kinv * blkdiag(Kxw{:})' * (CPhi_big * Y0); % S x length(nList)
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
        trKWW = trKWW + trace(Kinv{xIdx} * params.W(j).moment(currLatent(1):currLatent(2),currLatent(1):currLatent(2)));
    end
    lbterm = lbterm - (logdetK + trKWW);
    
    if inferX
        % Compute posterior covariance of X
        Kx = construct_Kx_smdlag(params, T);
        SigX = blkdiag(Kx{:}) - blkdiag(Kxw{:}) * blkdiag(Kinv{:}) ...
             * (blkdiag(Kw{:}) - SigW) * blkdiag(Kinv{:}) * blkdiag(Kxw{:})'; % (MpT x MpT), p-> M -> T
        SigX = 0.5 * (SigX + SigX');

        % Now, take only subsets of SigX that will have future use.
        % Restructure these portions to match the structure expected by other
        % functions.
        lat_block_idxs = get_block_idxs(repmat(numGroups*T,1,xDim));
        group_block_idxs = get_block_idxs(repmat(T,1,numGroups));

        % (xDim*numGroups) X (xDim*numGroups) Posterior covariance for each timepoint
        params.X.cov(j).Vsm = nan(mp,mp,T); % (M -> p) x (M -> p) x T
        for p1 = 1:xDim
            currLatent1 = lat_block_idxs{p1};
            for p2 = 1:xDim
                currLatent2 = lat_block_idxs{p2};
                SigX_j = SigX(currLatent1(1):currLatent1(2),currLatent2(1):currLatent2(2)); % (MT x MT), M -> T
                for m1 = 1:numGroups
                    currGroup1 = group_block_idxs{m1};
                    for m2 = 1:numGroups
                        currGroup2 = group_block_idxs{m2};
                        SigX_jm = SigX_j(currGroup1(1):currGroup1(2),currGroup2(1):currGroup2(2)); % (T x T)
                        params.X.cov(j).Vsm((m1-1)*xDim+p1,(m2-1)*xDim+p2,:) ...
                            = diag(SigX_jm);
                    end
                end
            end
        end

        % (numGroups*T) x (numGroups*T) Posterior covariance for each GP
        params.X.cov(j).VsmGP = nan(numGroups*T,numGroups*T,xDim); % (T -> M) x (T -> M) x p
        idx = (1:numGroups:numGroups*T);
        for xIdx = 1:xDim
            currLatent = lat_block_idxs{xIdx};
            SigX_j = SigX(currLatent(1):currLatent(2),currLatent(1):currLatent(2)); % (MT x MT), M -> T
            for m1 = 1:numGroups
                currGroup1 = group_block_idxs{m1};
                for m2 = 1:numGroups
                    currGroup2 = group_block_idxs{m2};
                    SigX_jm = SigX_j(currGroup1(1):currGroup1(2),currGroup2(1):currGroup2(2)); % (T x T)
                    params.X.cov(j).VsmGP(idx+(m1-1),idx+(m2-1),xIdx) ...
                        = SigX_jm;
                end
            end
        end

        % Compute posterior mean of X
        xsmMat = blkdiag(Kxw{:}) * blkdiag(Kinv{:}) * wsmMat; % (MpT x length(nList)), p -> M -> T
        % Restructure X to match the structure expected in other functions
        for i = 1:length(nList)
            n = nList(i);
            seq(n).xsm = zeros(numGroups*xDim,T);
            for xIdx = 1:xDim
                currLatent = lat_block_idxs{xIdx};
                xsm_j = xsmMat(currLatent(1):currLatent(2),i); % (MT x 1), M -> T
                for groupIdx = 1:numGroups
                    currGroup = group_block_idxs{groupIdx};
                    xsm_jm = xsm_j(currGroup(1):currGroup(2)); % (T x 1)
                    seq(n).xsm((groupIdx-1)*xDim+xIdx,:) = xsm_jm'; % (Mp x 1), M -> p
                end
            end
        end
    end
end
