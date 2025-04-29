function precomp = makePrecomp_GPdelays_freq(seq,params,CPhiC,Y0)
%
% precomp = makePrecomp_GPdelays_freq(seq,params,CPhiC,Y0)
%
% Description: Precompute posterior second moments for the GP delay
%              update (frequency domain approximation).
%
% Arguments:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT of the observed data
%                     xfft (xDim x T) -- unitary FFT of the latent 
%                                        posterior mean at each frequency
%
%     params -- Structure containing mDLAG model parameters.
%               Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         X.spec     -- data structure whose jth entry,
%                       corresponding to a group of trials of the 
%                       same length, has fields
%                           T       -- int; number of time steps for
%                                      this trial group
%                           Sx_post -- (xDim x xDim x T) array; 
%                                      posterior spectrum at each 
%                                      frequency
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
%
%     CPhiC  -- (1 x numGroups) cell array; CPhiC{m} is a (xDim x xDim)
%               array containing the expectation of C_m'*Phi_m*C_m
%     Y0     -- (1 x length(Tu)) cell array; Y0{j} contains the 
%               mean-subtracted observations for all trials of length Tu(j)
%
% Outputs
%     precomp -- Structure whose m-th entry contains precomputations for 
%                the m-th observation group:
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        T -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        XXCPhiC -- (xDim x xDim x T) array; An 
%                                   intermediate quantity to be reused.
%                        YXPhiC -- (yDims(m) x xDim x T) array; a common 
%                                  intermediate quantity to be reused  
%                NOTE: The first entry of precomp is empty because 
%                      delays to the first group are not updated.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Aug 2023 -- Initial full revision. 
%     29 Aug 2024 -- Added Y0 as an argument.

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
obs_block_idxs = get_block_idxs(yDims);

% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);

% Initialize output structure
for j = 1:length(Tu)
    for groupIdx = 1:numGroups
        T = Tu(j);
        nList = find(Tall == T);
        precomp(groupIdx).Tu(j).T = T;
        precomp(groupIdx).Tu(j).numTrials = length(nList);
        if groupIdx <= 1
            precomp(groupIdx).Tu(j).XXCPhiC = [];
            precomp(groupIdx).Tu(j).YXPhiC = [];
        else
            precomp(groupIdx).Tu(j).XXCPhiC = zeros(xDim,xDim,T);
            precomp(groupIdx).Tu(j).YXPhiC = zeros(yDims(groupIdx),xDim,T);
        end
    end
end

% Fill in output structure
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    X = cat(3,seq(nList).xfft); % (xDim x T x N)

    % XX
    XX = length(nList).*params.X.spec(j).Sx_post ...
        + sum(repmat(permute(X,[1 4 2 3]),1,xDim,1,1) ...
        .* repmat(permute(conj(X),[4 1 2 3]),xDim,1,1,1),4); % (xDim x xDim x T)

    for groupIdx = 2:numGroups
        currObsGroup = obs_block_idxs{groupIdx};
        obsIdxs = currObsGroup(1):currObsGroup(2); % Current observation group

        % XXCPhiC
        precomp(groupIdx).Tu(j).XXCPhiC = permute(XX,[2 1 3]) .* repmat(CPhiC{groupIdx},1,1,T);

        % XY
        XY = sum(repmat(permute(X,[1 4 2 3]),1,yDims(groupIdx),1,1) ...
             .* repmat(permute(conj(Y0{j}(obsIdxs,:,:)),[4 1 2 3]),xDim,1,1,1),4); % (xDim x yDims(groupIdx) x T)

        % YXPhiC
        PhiC = params.phi.mean(obsIdxs) .* params.C.means{groupIdx};
        precomp(groupIdx).Tu(j).YXPhiC = permute(XY,[2 1 3]) .* repmat(PhiC,1,1,T);
    end
end
