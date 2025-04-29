
function precomp = makePrecomp_GPtimescales_freq(seq,params)
%
% precomp = makePrecomp_GPtimescales_freq(seq,params)
%
% Description: Precompute posterior second moment for the GP timescale
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
% Outputs
%     precomp -- Structure whose ith entry contains precomputations for 
%                the i-th latent state:
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        T -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        XX -- (T x 1) array; Posterior second moment
%                              of this latent state        
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Aug 2023 -- Initial full revision.
%     07 Aug 2024 -- Only diagonal elements of XX are needed, so lower
%                    bound and gradients can be computed much more
%                    efficiently.

% Initialize relevant variables
xDim = params.xDim;

% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);

% Initialize output structure
for j = 1:length(Tu)
    for k = 1:xDim
        T = Tu(j);
        nList = find(Tall == T);
        precomp(k).Tu(j).T = T;
        precomp(k).Tu(j).numTrials = length(nList);
        precomp(k).Tu(j).XX = zeros(T, 1);
    end
end

% Fill in output structure
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    X = permute(cat(3,seq(nList).xfft), [2 3 1]);  % (T x N x xDim)
    for k = 1:xDim
        % XX
        Sx_post = reshape(params.X.spec(j).Sx_post(k,k,:),[],1);
        precomp(k).Tu(j).XX ...
            = precomp(k).Tu(j).numTrials.*Sx_post ...
            + sum(X(:,:,k) .* conj(X(:,:,k)), 2);
    end
end
