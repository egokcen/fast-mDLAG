function seq = generate_latents_mdlag_freq(params, T, N, varargin)
% 
% seq = generate_latents_mdlag_freq(params, T, N,...)
%
% Description: Generate N independent sequences of length T samples, 
%              according to a zero-mean Gaussian process defined by the 
%              mDLAG state model. Sequences are generated approximately in 
%              the frequency domain. The true covariances are thus 
%              circulant. Only the complex-conjugate portion of the
%              frequency domain latents should be considered.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance. The 
%                               following GP kernels are currently 
%                               supported:
%                                   'rbf'    -- Radial basis function, or 
%                                               squared exponential
%                                   'sg'     -- Spectral Gaussian
%                                   'exp'    -- Exponential
%                                   'expcos' -- Exponential-cosine
%                    The following fields depend on GP covariance type:
%                        For 'rbf':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma)                                                    
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'sg':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma) 
%                            nu    -- (1 x xDim) array; center frequencies; 
%                                     convert to 1/time via nu./binWidth 
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'exp':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ gamma                                                    
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'expcos':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ gamma 
%                            nu    -- (1 x xDim) array; center frequencies; 
%                                     convert to 1/time via nu./binWidth 
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                    Cs      -- (1 x numGroups) cell array; List of factor 
%                               loadings 
%                               {(y1Dim x xDim), (y2Dim x xDim), ...}
%                    alphas  -- (numGroups x xDim) array; alpha parameter 
%                               values
%                    phis    -- (1 x numGroups) cell array; List of 
%                               observation precision parameters 
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    ds      -- (1 x numGroups) cell array; List of data 
%                               means
%                               {(y1Dim x 1), (y2Dim x 1), ...}
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%     T         -- (1 x N) int array; T(n) gives the number of samples 
%                  (time length) for sequence n. If T is a scalar
%                  (length-1), then all sequences will be length-T.
%     N         -- int; number of sequences
%
%     Optional:
%
%     freqfield -- string; field name for frequency domain latents
%                  (default: 'xfft')
%     timefield -- string; field name for (time-delayed) time domain 
%                  latents (default: 'xsm')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (freqfield) -- (xDim x T) array; frequency domain latents
%                (timefield) -- (numGroups*xDim x T) array; (time-delayed)
%                               time domain latents
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     03 Sep 2023 -- Initial full revision.
%     07 Mar 2024 -- Updated documentation to reflect new GP kernels.
%     13 Jun 2024 -- Updated scaling in frequency domain so that time
%                    domain signals have unit variance.
%     12 Aug 2024 -- Updated with make_S_mdlag, which now returns only
%                    a vector of diagonal elements.

freqfield = 'xfft';
timefield = 'xsm';
assignopts(who, varargin);
   
% If T is a scalar, then give all sequences the same length
if length(T) <= 1
    T = repmat(T,1,N);
end

% Group trials of same length together
Tu = unique(T);

% Initialize output structure
xDim = params.xDim;
for n = 1:N
    seq(n).trialId = n;
    seq(n).T = T(n);
    seq(n).(freqfield) = nan(xDim,T(n));
end
    
% Generate all trials of the same length
for j = 1:length(Tu)
    Tj = Tu(j);
    nList = find(T == Tj);
    N = length(nList);

    % Handle even and odd sequence lengths
    freqs = (-floor(Tj/2):floor((Tj-1)/2))./Tj;

    X = zeros(xDim,Tj,N);
    for f = 1:Tj
        % Construct S, the spectral density matrix
        Sf = make_S_mdlag(params,freqs(f));
        % Generate latents in the frequency domain
        % NOTE: Normally, a complex-valued normal random variable with 
        %       variance sigma would be generated according to
        %           sqrt(sigma/2) * (randn + 1i*randn)
        %       Here, we're generating a complex-valued signal in the
        %       frequency domain, and later we'll take only the real part
        %       of its time domain counterpart. Equivalently, we'll be 
        %       taking only the conjugate-symmetric portion of the
        %       frequency-domain signal. That step will throw out half of
        %       the signal's total power. We want signals in the time 
        %       domain, however, to have unit variance. Therefore, we'll 
        %       effectively double the variance in this next step, and
        %       remove the factor of sqrt(1/2). This choice will give us
        %       unit variance in the time domain.
        for k = 1:xDim
            X(k,f,:) = sqrt(Sf(k)).*(randn(1,1,N) + 1i.*randn(1,1,N));
        end
    end
    % Collect latents into output structure
    for n = 1:N
        seq(nList(n)).(freqfield) = X(:,:,n);
    end
end

% Convert latents to the time domain
seq = freq2time_mdlag(seq, params, ...
                      'infield', freqfield, ...
                      'outfield', timefield);
