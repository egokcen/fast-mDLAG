function [seq, params, lbterms] = inferX_freq(seq, params, varargin)
%
% [seq, params, lbterms] = inferX_freq(seq, params, ...)
%
% Description: Infer latent variables X given observations Y (in seq)
%              and mDLAG model parameters in params, using an approximate 
%              frequency domain approach.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId         -- unique trial identifier
%                    T (1 x 1)       -- number of timesteps
%                    yfft (yDim x T) -- unitary FFT of observed data
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
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
%    Optional:
%   
%    CPhiC        -- (1 x numGroups) cell array; CPhiC{m} is a 
%                    (xDim x xDim) array containing the precomputed 
%                    expectation of C_m'*Phi_m*C_m (default: {}, in which 
%                    case CPhiC will be computed from scratch)
%    Y0           -- (1 x length(Tu)) cell array; Y0{j} contains the 
%                    mean-subtracted observations for all trials of
%                    length Tu(j) (default: {}, in which case Y0 will be 
%                    computed from scratch)
%
% Outputs:
%
%     seq    -- data structure whose nth entry has new fields (these fields 
%               are added to existing fields in the seq input argument)
%                 xfft  -- (xDim x T) array; unitary FFT of the posterior 
%                           mean at each timepoint
%     params -- Structure containing mDLAG model parameters with new field
%                 X.spec -- data structure whose jth entry,
%                              corresponding to a group of trials of the 
%                              same length, has fields
%                               T       -- int; number of time steps for
%                                          this trial group
%                               Sx_post -- (xDim x xDim x T) array; 
%                                          posterior spectrum at each 
%                                          frequency
%                             NOTE: For mDLAG, posterior 
%                                   covariances/spectra of X are the same 
%                                   for trials of the same length.
%     lbterms -- Structure containing terms that contribute to the
%                variational lower bound:
%                  logdet_Sx_post  -- float; log|Sx_post|
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Aug 2023 -- Initial full revision.
%     12 Aug 2024 -- Made several runtime optimizations.
%     29 Aug 2024 -- Added option to compute Y0 or pass it as an argument.

CPhiC = {};
Y0 = {};
assignopts(who, varargin);

% Initialize other relevant variables
yDims = params.yDims; 
yDim = sum(yDims);
xDim = params.xDim;
numGroups = length(yDims);
block_idxs = get_block_idxs(yDims);

% Compute C'*Phi, an intermediate term
CPhi = cell(1,numGroups);
for groupIdx = 1:numGroups
    currObsGroup = block_idxs{groupIdx};
    obsIdxs = currObsGroup(1):currObsGroup(2); % Current observation group
    CPhi{groupIdx} = params.C.means{groupIdx}.' .* params.phi.mean(obsIdxs).';
end

if isempty(CPhiC)
    CPhiC = computeCPhiC(params);
end

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

if isempty(Y0)
    Y0 = cell(1,length(Tu));
end

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.
lbterms.logdet_Sx_post = 0;
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);

    % Handle even and odd sequence lengths
    freqs = (-floor(T/2):floor((T-1)/2))./T;

    if isempty(Y0{j})
        % Find the index of zero frequency
        zeroIdx = floor(T/2)+1;
        % Zero-center the data
        Y0{j} = cat(3,seq(nList).yfft);
        Y0{j}(:,zeroIdx,:) = Y0{j}(:,zeroIdx,:) ...
            - repmat(sqrt(T).*params.d.mean,[1,1,length(nList)]);
    end

    % Initialize output structures
    xfftMat = zeros(xDim,length(nList),T);
    params.X.spec(j).T = T;
    params.X.spec(j).Sx_post = nan(xDim,xDim,T);
    for n = nList
        seq(n).xfft = nan(xDim,T); 
    end
    
    for f = 1:T
        % Construct S, the spectral density matrix
        Sx = make_S_mdlag(params,freqs(f));  % (xDim x 1)
        Sx_inv = 1./Sx;  % (xDim x 1)

        % Construct Q, the time-delay operator
        Q = exp(-1i*2*pi*freqs(f).*params.D);  % numGroups x xDim
        
        % Posterior spectral density matrix
        Sx_post_inv = diag(Sx_inv);  % xDim x xDim
        for groupIdx = 1:numGroups
            Sx_post_inv = Sx_post_inv + Q(groupIdx,:)' .* CPhiC{groupIdx} .* Q(groupIdx,:);
        end
        Sx_post_inv = (Sx_post_inv + Sx_post_inv')./2; % Ensure Hermitian
        Sx_post = inv(Sx_post_inv);   % xDim x xDim
        params.X.spec(j).Sx_post(:,:,f) = Sx_post;
        
        % Posterior mean of X
        for groupIdx = 1:numGroups
            currObsGroup = block_idxs{groupIdx};
            xfftMat(:,:,f) = xfftMat(:,:,f) + (Q(groupIdx,:)' .* CPhi{groupIdx}) ...
                * permute(Y0{j}(currObsGroup(1):currObsGroup(2),f,:), [1 3 2]);
        end
        xfftMat(:,:,f)  = Sx_post * xfftMat(:,:,f);  % xDim x length(nList)
        
        % Lower bound term
        lbterms.logdet_Sx_post = lbterms.logdet_Sx_post + length(nList)*logdet(Sx_post);

    end

    % Reshape X and collect it in seq
    ctr = 1;
    for n = nList
      seq(n).xfft(:,:) = permute(xfftMat(:,ctr,:), [1 3 2]);      
      ctr = ctr + 1;
    end

end
