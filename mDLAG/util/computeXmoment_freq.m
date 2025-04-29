function [Xmoment, Xdelayed] = computeXmoment_freq(seq,params)
%
% [Xmoment, Xdelayed] = computeXmoment_freq(seq,params)
%
% Description: Helper function to compute the (time-delayed) posterior 
%              second moments (in the approximate frequency domain) of X 
%              for each group.
%
% Arguments:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT of observed data
%                     xfft (xDim x T) -- unitary FFT of the posterior mean
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
%         X.spec -- data structure whose jth entry,
%                   corresponding to a group of trials of the 
%                   same length, has fields
%                    T       -- int; number of time steps for
%                               this trial group
%                    Sx_post -- (xDim x xDim x T) array; 
%                               posterior spectrum at each 
%                               frequency
%                  NOTE: For mDLAG, posterior 
%                        covariances/spectra of X are the same 
%                        for trials of the same length.
%         covType -- string; type of GP covariance. The 
%                    following GP kernels are currently 
%                    supported:
%                        'rbf'    -- Radial basis function, or 
%                                    squared exponential
%                        'sg'     -- Spectral Gaussian
%                        'exp'    -- Exponential
%                        'expcos' -- Exponential-cosine
%         The following fields depend on GP covariance type:
%             For 'rbf':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)                                                      
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'sg':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)    
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'exp':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma                                                      
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'expcos':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma    
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
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
%
%     Xmoment  -- (1 x numGroups) cell array; Xmoment{m} is a (xDim x xDim)
%                 array containing second moments of X for the mth group.  
%     Xdelayed -- data structure whose jth entry, corresponding to a 
%                 group of trials of the same length, has fields
%                     T  -- int; number of time steps for this trial group
%                     QX -- (1 x numGroups) cell array; QX{m} is a 
%                           (xDim x T x N) array containing the 
%                           (time-delayed) posterior mean of X across
%                           trials in this trial group
%                 Xdelayed is reused elsewhere in the EM algorithm.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Aug 2023 -- Initial full revision. 
%     12 Aug 2024 -- Made several runtime optimizations.

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
Tall = [seq.T];

% Find unique numbers of trial lengths
Tu = unique(Tall);

% Initialize Xmoment
Xmoment = cell(1,numGroups);
for groupIdx = 1:numGroups
    Xmoment{groupIdx} = zeros(xDim);
end
% Loop over each trial length (each of Tu)
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    N = length(nList);
    X = cat(3,seq(nList).xfft); % (xDim x T x N)

    % Handle even and odd sequence lengths
    freqs = (-floor(T/2):floor((T-1)/2))./T;

    % Initialize Xdelayed for the current trial group
    Xdelayed(j).T = T;
    Xdelayed(j).QX = cell(1,numGroups);

    for groupIdx = 1:numGroups
        % Construct Qm, the time-delay operator for the current group
        Qm = exp(-1i*2*pi*params.D(groupIdx,:)'*freqs);
        % Compute time-delayed latent means
        QX = repmat(Qm,1,N) .* reshape(X,xDim,[]);        % (xDim x T*N)
        Xdelayed(j).QX{groupIdx} = reshape(QX,xDim,T,N);  % (xDim x T x N)
        % Second moment of X
        Qm_big = repmat(permute(Qm,[1 3 2]),[1 xDim 1]);
        Qm_herm_big = permute(conj(Qm_big),[2 1 3]);
        Sm = N.*(sum(Qm_big.*params.X.spec(j).Sx_post.*Qm_herm_big,3));
        Xmoment{groupIdx} = Xmoment{groupIdx} + Sm + QX * QX';
    end
end
