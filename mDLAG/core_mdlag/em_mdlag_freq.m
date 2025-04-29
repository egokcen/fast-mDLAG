function [currentParams,seq,trackedParams,flags] ...
    = em_mdlag_freq(currentParams,seq,xDim,varargin)
%
% [currentParams,seq,trackedParams,flags] = em_mdlag_freq(currentParams,seq,xDim, ...)
%
% Description: Fit a (multi-group) Delayed Latents Across Groups (mDLAG) 
%              model using a variational EM algorithm with mean-field
%              approximation. For additional speedup, use an approximate
%              frequency domain approach.
%
% Arguments:
%
%     Required:
%
%     currentParams  -- Structure containing mDLAG model parameters.
%                       Contains the fields
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
%     seq     -- data structure, whose nth entry (corresponding to the nth 
%                trial) has fields
%                    trialId         -- unique trial identifier
%                    T (1 x 1)       -- number of timesteps
%                    y (yDim x T)    -- observed data
%                    yfft (yDim x T) -- unitary FFT of observed data
%     xDim    -- int; number of latent variables
%
%     Optional:
%
%     prior     -- structure with the following fields:
%                    d.beta  -- positive float; precision of mean parameter
%                               generative model (Gaussian)
%                    phi.a   -- positive float; 'a' shape parameter of
%                               observation precision (phi) generative 
%                               model (Gamma with mean a/b)
%                    phi.b   -- positive float; 'b' scale parameter of
%                               observation precision (phi) generative
%                               model (Gamma with mean a/b)
%                    alpha.a -- positive float; 'a' shape parameter of 
%                               alpha parameter generative model 
%                               (Gamma with mean a/b)
%                    alpha.b -- positive float; 'b' scale parameter of
%                               alpha parameter generative model 
%                               (Gamma with mean a/b)
%                    (default: 1e-12 for all values)
%
%     tol       -- float; stopping criterion for EM (default: 1e-8)
%     maxIters  -- int; maximum number of EM iterations (default: 1e6)
%     freqLB    -- int; lower bound is computed every freqLB EM iterations.
%                  freqLB = 1 means that the lower bound is computed every
%                  iteration. (default: 1)
%     freqParam -- int; store intermediate values for delays and timescales
%                  every freqParam EM iterations (default: 10)
%     learnDelays -- logical; set true to learn delay parameters;
%                    otherwise, delays will remain fixed at their initial
%                    value. (default: true)
%     learnTimescales -- logical; set true to learn timescale parameters;
%                        otherwise, timescales will remain fixed at their 
%                        initial value. (default: true)
%                        NOTE: Only 'rbf' and 'exp' kernels are currently
%                              supported.
%     verbose   -- boolean; specifies whether to display status messages
%                  (default: false)
%     minVarFrac -- float; fraction of overall data variance for each 
%                   observed dimension to set as the private variance 
%                   floor. (default: 1e-3)
%     maxDelayFrac -- float in range [0,1]; Constrain estimated delays to
%                     be no more than a certain fraction of the trial 
%                     length. (default: 0.5)
%     maxTauFrac -- float in range [0,1]; Constrain estimated timescales to
%                   be no more than a certain fraction of the trial 
%                   length. (default: 1.0)
%     trackedParams -- structure containing the tracked parameters from a
%                      previous fitting attempt, with the intention of 
%                      starting where that attempt left off. See 
%                      trackedParams below. (default: {})
%     pruneX    -- logical; set true to remove dimensions from X that
%                  become inactive in all groups. Can speed up EM runtime
%                  and improve memory efficiency. (default: true)
%     saveXcov  -- logical; set true to save posterior covariance of latent
%                  variables X. For large datasets, these structures can
%                  use a lot of memory. (default: false)
%     saveCcov  -- logical; set true to save posterior covariance and 
%                  of C. For large yDim and xDim, these structures can use
%                  a lot of memory. (default: false)
%
% Outputs:
%
%     currentParams -- Structure containing mDLAG model parameters returned 
%                  by EM algorithm (same format as currentParams) with the
%                  addition of (if saveX is true)
%                    X.spec -- data structure whose jth entry,
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
%     seq       -- data structure whose nth entry has new fields 
%                  (these fields are added to existing fields in the seq 
%                  input argument)
%                  xfft  -- (xDim x T) array; unitary FFT of the posterior 
%                           mean at each timepoint
%     trackedParams -- structure containing parameters tracked throughout 
%                      fitting:
%       lb        -- (1 x numIters) array; variational lower bound at each 
%                    iteration
%       iterTime  -- (1 x numIters) array; computation time for each EM
%                    iteration
%       alphas   -- (1 x numIters) cell arry; estimated ARD parameters
%                   (alpha.mean) after each EM iteration.
%       gp_params -- structure tracking kernel-dependent GP parameters. 
%                    Contains the following fields:
%         For 'rbf' or 'exp':
%           Ds       -- (1 x numIters) cell array; the estimated delay 
%                       matrix (D) after each EM iteration.
%           gams     -- (1 x numIters) cell arry; estimated gamma after 
%                       each EM iteration.
%         For 'sg' or 'expcos':
%           Ds       -- (1 x numIters) cell array; the estimated delay 
%                       matrix (D) after each EM iteration.
%           gams     -- (1 x numIters) cell arry; estimated gamma after 
%                       each EM iteration.
%           nus      -- (1 x numIters) cell arry; estimated nu after
%                       each EM iteration.
%     flags     -- structure containing various status messages about
%                  the fitting procedure:
%         Convergence          -- logical; 1 if lower bound converged 
%                                 before reaching maxIters EM iterations; 
%                                 0 otherwise
%         DecreasingLowerBound -- logical; 1 if lower bound decreased 
%                                 during fitting; 0 otherwise.
%         PrivateVarianceFloor -- logical; 1 if private variance floor was 
%                                 used on any values of phi; 0 otherwise.
%         xDimsRemoved         -- int; Number of latent dimensions
%                                 removed (if pruneX is true) due to 
%                                 low variance in all groups.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Aug 2023 -- Initial full revision.
%     03 Sep 2023 -- Added ability to handle alternative GP kernels.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.
%     18 Jun 2024 -- Added learnTimescales option for rbf, exp kernels.
%     12 Aug 2024 -- Revised with a few runtime optimizations.
%     29 Aug 2024 -- Making better reuse of Y0 data structure.

prior.d.beta = 1e-12;
prior.phi.a = 1e-12;
prior.phi.b = 1e-12;
prior.alpha.a = 1e-12;
prior.alpha.b = 1e-12;
tol = 1e-8; 
maxIters = 1e6;
freqLB = 1;
freqParam = 10;
learnDelays = true;
learnTimescales = true;
verbose = false;
minVarFrac = 1e-3;
maxDelayFrac  = 0.5;
maxTauFrac = 1.0;
trackedParams = {};
pruneX = true;
saveXcov = false;
saveCcov = false;
assignopts(who, varargin);

% Initialize common variables
yDims = currentParams.yDims;
numGroups = length(yDims);
% Useful for extracting correct-sized blocks from matrices later
block_idxs = get_block_idxs(yDims);

% If currentParams was initialized to include C.covs, but C.covs should
% not be saved, then remove it.
if ~saveCcov && isfield(currentParams.C, 'covs')
    currentParams.C = rmfield(currentParams.C, 'covs');
end

% If C.covs should be saved, but currentParams was not initialized with it,
% then initialize the structure (initial value doesn't matter here).
if saveCcov && ~isfield(currentParams.C, 'covs')
    currentParams.C.covs = currentParams.C.moments; 
end

% Get features of the joint data matrix
[yDim, NT] = size([seq.yfft]);
sum_Y2 = sum(abs([seq.yfft]).^2,2);  % To be re-used later

% Get private variance floors
varFloor = minVarFrac * var([seq.y],0,2);

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

% Group the data by trials of the same length
Y = cell(1,length(Tu));
Y0 = cell(1,length(Tu));  % Mean-centered observations. Initialize for now.
sum_Yf0 = cell(1,length(Tu));  % Sum over trials of zero-frequency data. 
for j = 1:length(Tu)
    % Process all trials with length T
    T = Tu(j);
    nList = find(Tall == T);
    Y{j} = cat(3,seq(nList).yfft);        % (yDim x T x N)
    % Find the index of zero frequency
    zeroIdx = floor(T/2)+1;
    % Zero-centered observations
    Y0{j} = Y{j};
    Y0{j}(:,zeroIdx,:) = Y0{j}(:,zeroIdx,:) ...
        - repmat(sqrt(T).*currentParams.d.mean,[1,1,length(nList)]);
    % Zero-frequency observations
    sum_Yf0{j} = sum(Y{j}(:,zeroIdx,:), 3);  % (yDim x 1)
end

% Define constants
% Latent magnitude must be greater than this value to remain in the model
PRUNEX_TOL = 1e-6;

% Constant factors in elbo
const_lik = -((yDim*NT)/2) * log(2*pi);

% Constants in observation noise model
alogb_phi = prior.phi.a * log(prior.phi.b);
loggammaa_phi_prior = gammaln(prior.phi.a);
loggammaa_phi_post = gammaln(currentParams.phi.a);
digammaa_phi = psi(currentParams.phi.a);

% Constants in ARD model
alogb_alpha = prior.alpha.a * log(prior.alpha.b);
loggammaa_alpha_prior = gammaln(prior.alpha.a);
loggammaa_alpha_post = gammaln(currentParams.alpha.a); % (numGroups x 1) array
digammaa_alpha = psi(currentParams.alpha.a);           % (numGroups x 1) array 

% Make sure initial delays are within specified constraints
% Convert maxDelayFrac to units of "time steps", the same units as in D
constraints.maxDelay = maxDelayFrac*min([seq.T]); 

% If delays are outside the range (minDelay,maxDelay), then clamp them
% within the range
currentParams.D(currentParams.D >= constraints.maxDelay) = 0.95*constraints.maxDelay;
currentParams.D(currentParams.D <= -constraints.maxDelay) = -0.95*constraints.maxDelay;

% Convert maxTauFrac to unitless quantity 'gamma'
switch currentParams.covType
    case {'rbf', 'sg'}
        constraints.minGamma = 1/(maxTauFrac*min([seq.T]))^2;
    case {'exp', 'expcos'}
        constraints.minGamma = 1/(maxTauFrac*min([seq.T]));
end

if isempty(trackedParams)
    % Start convergence tracking fresh
    lbi = -Inf;                   % Initial variational lower bound
    lb = [];                      % Lower bound at each iteration
    iterTime = [];                % Time it takes to complete each iteration
    gp_params.Ds = {currentParams.D}; % Delays each iteration
    switch currentParams.covType
        case {'rbf', 'exp'}     
            gp_params.gams = {currentParams.gamma}; % Timescales each iteration
        case {'sg', 'expcos'}
            gp_params.gams = {currentParams.gamma}; % Timescales each iteration
            gp_params.nus = {currentParams.nu}; % Center frequencies each iteration
    end
    alphas = {currentParams.alpha.mean}; % ARD parameters each iteration
    startIter = 1;                % Initial value for EM loop index
else
    % Start convergence tracking based on where a previous attempt left
    % off.
    lb             = trackedParams.lb;       % Lower bound at each iteration
    lbi            = trackedParams.lb(end);  % Most recent lower bound
    lbbase         = trackedParams.lb(2);    % Base lower bound
    iterTime       = trackedParams.iterTime; % Time it takes to complete each iteration
    gp_params.Ds   = trackedParams.gp_params.Ds; % Delays each iteration
    switch currentParams.covType
        case {'rbf', 'exp'}
            gp_params.gams = trackedParams.gp_params.gams; % Timescales each iteration
        case {'sg', 'expcos'}
            gp_params.gams = trackedParams.gp_params.gams; % Timescales each iteration
            gp_params.nus  = trackedParams.gp_params.nus;  % Center frequencies each iteration
    end
    alphas    = trackedParams.alphas;   % ARD parameters each iteration
    startIter = length(lb) + 1;         % Initial value for EM loop index
end

% Initialize status flags
flags.Convergence = 0;
flags.DecreasingLowerBound = 0;
flags.PrivateVarianceFloor = 0;
flags.xDimsRemoved = 0;

% Initialize posterior mean of latent states X
[seq,~,~] = inferX_freq(seq, currentParams);

% Begin EM iterations
for i = startIter:maxIters
    
    % Determine when to actually compute the lower bound
    if (rem(i, freqLB) == 0) || (i<=2) || (i == maxIters)
        getLB = true;
    else
        getLB = false;
    end
    
    tic; % For tracking the computation time of each iteration
    
    % ======================
    % Latent states, X
    % ======================
    
    % Check if any latent states need to be removed
    if pruneX
        % Concatenate latent states across time and trials
        X = [seq.xfft];
        % To be kept, latent states must show at least minimal activity
        keptXdims = find(mean(abs(X).^2,2) > PRUNEX_TOL);
        if length(keptXdims) < xDim
            currentParams = getSubsetXDims_params(currentParams,keptXdims);
            [gp_params,alphas] ...
                = getSubsetXDims_trackedParams(currentParams.covType,...
                                               gp_params,...
                                               alphas,...
                                               keptXdims);
            flags.xDimsRemoved = flags.xDimsRemoved + (xDim - currentParams.xDim);
            xDim = currentParams.xDim;
            if xDim <= 0
                % Stop fitting if no significant latent dimensions remain.
                break;
            end
        end
    end

    % Precompute the expectation of C'*Phi*C, an intermediate term
    % shared across X and GP parameter updates
    CPhiC = computeCPhiC(currentParams);
    
    % Posterior mean and portions of the posterior covariance of X
    [seq, currentParams, lbterms] = inferX_freq(seq, currentParams, 'CPhiC', CPhiC, 'Y0', Y0);
    
    % ======================
    % GP parameters
    % ======================
    
    res = learnGPparams_mdlag_freq(seq, currentParams, constraints, CPhiC, Y0, ...
                                   'learnDelays', learnDelays, ...
                                   'learnTimescales', learnTimescales);
    switch currentParams.covType
        case {'rbf', 'exp'}
            currentParams.gamma = res.gamma; 
            if numGroups > 1
                currentParams.D = res.D;
            end
            if (rem(i, freqParam) == 0) || (i == maxIters)
                % Store current delays and timescales and compute
                % change since last computation
                gp_params.Ds = [gp_params.Ds {currentParams.D}];
                gp_params.gams = [gp_params.gams {currentParams.gamma}];
            end
        case {'sg', 'expcos'}
            currentParams.gamma = res.gamma; 
            currentParams.nu = res.nu;
            if numGroups > 1
                currentParams.D = res.D;
            end
            if (rem(i, freqParam) == 0) || (i == maxIters)
                % Store current delays and timescales and compute
                % change since last computation
                gp_params.Ds = [gp_params.Ds {currentParams.D}];
                gp_params.gams = [gp_params.gams {currentParams.gamma}];
                gp_params.nus = [gp_params.nus {currentParams.nu}];
            end
    end

    % (Time-delayed) posterior first and second moments of X for each group
    [XX, Xdelayed] = computeXmoment_freq(seq,currentParams);
    
    % ======================
    % Mean parameter, d
    % ======================
    
    currentParams.d.cov = 1./(prior.d.beta + NT*currentParams.phi.mean);
    currentParams.d.mean = 0;
    for j = 1:length(Tu)
        % Process all trials with length T
        T = Tu(j);
        nList = find(Tall == T);
        % Find the index of zero frequency
        zeroIdx = floor(T/2)+1;
        % We only need the zero-frequency components for this computation
        X = cat(3,seq(nList).xfft);            % (xDim x T x N)
        X = permute(X(:,zeroIdx,:), [1 3 2]);  % (xDim x N)
        currentParams.d.mean = currentParams.d.mean ...
            + sqrt(T).*sum((permute(Y{j}(:,zeroIdx,:), [1 3 2]) ...
            - blkdiag(currentParams.C.means{:}) * repmat(X,numGroups,1)),2);
    end
    currentParams.d.mean = currentParams.d.cov .* currentParams.phi.mean ...
        .* currentParams.d.mean;
    dd = diag(diag(currentParams.d.cov) ...
        + currentParams.d.mean*currentParams.d.mean'); % Second moment of d

    % Update zero-centered observations
    for j = 1:length(Tu)
        % Process all trials with length T
        T = Tu(j);
        % Find the index of zero frequency
        zeroIdx = floor(T/2)+1;
        % Zero-center the data
        Y0{j}(:,zeroIdx,:) = Y{j}(:,zeroIdx,:) - sqrt(T).*currentParams.d.mean;
    end
    
    % ======================
    % Loading matrices, C
    % ======================
    logdetC = 0; % To be used in calculation of lower bound
    XY = cell(1,numGroups); % Inner product to be reused in phi update
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_m = currentParams.phi.mean(currGroup(1):currGroup(2));
        alpha_m = diag(currentParams.alpha.mean(groupIdx,:));  % (xDim x xDim) array
        XY{groupIdx} = zeros(xDim, yDims(groupIdx));
        for j = 1:length(Tu)
            % Process all trials with length T
            XY{groupIdx} = XY{groupIdx} + reshape(Xdelayed(j).QX{groupIdx},xDim,[]) ...
                * reshape(Y0{j}(currGroup(1):currGroup(2),:,:),yDims(groupIdx),[])';
        end
        for yIdx = 1:yDims(groupIdx)
            % Covariance
            covC = inv(alpha_m + phi_m(yIdx) .* real(XX{groupIdx}));
            covC = 0.5 * (covC + covC'); % Ensure symmetry
            logdetC = logdetC + logdet(covC);
            if saveCcov
                currentParams.C.covs{groupIdx}{yIdx} = covC;
            end
            % Mean
            currentParams.C.means{groupIdx}(yIdx,:) ...
                = (phi_m(yIdx) * covC * real(XY{groupIdx}(:,yIdx)))';
            % Second moment
            currentParams.C.moments{groupIdx}{yIdx} = covC ...
                + currentParams.C.means{groupIdx}(yIdx,:)' * currentParams.C.means{groupIdx}(yIdx,:);
        end
    end
    
    % ======================
    % ARD parameters, alpha
    % ======================
    
    for groupIdx = 1:numGroups
        % Second moments of C_m (m = groupIdx)
        CC_m = sum(cat(3,currentParams.C.moments{groupIdx}{:}),3);
        currentParams.alpha.b(groupIdx,:) = (prior.alpha.b + diag(CC_m)./2)';
        currentParams.alpha.mean(groupIdx,:) ...
            = currentParams.alpha.a(groupIdx) ./ currentParams.alpha.b(groupIdx,:); 
    end
    if (rem(i, freqParam) == 0) || (i == maxIters)
        alphas = [alphas {currentParams.alpha.mean}];
    end
    
    % ======================
    % Noise precisions, phi
    % ======================
    
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        phi_bm = repmat(prior.phi.b,yDims(groupIdx),1);
        dd_m = dd(currGroup(1):currGroup(2));
        d_m = currentParams.d.mean(currGroup(1):currGroup(2));
        sum_Y2_m = sum_Y2(currGroup(1):currGroup(2));

        % Handle zero-frequency terms
        for j = 1:length(Tu)
            % Process all trials with length T
            T = Tu(j);
            nList = find(Tall == T);
            N = length(nList);
            phi_bm = phi_bm + 0.5 * real(N*T * dd_m ...
                - 2*sqrt(T).*sum_Yf0{j}(currGroup(1):currGroup(2)) .* d_m);
        end

        % Handle the remaining terms
        phi_bm = phi_bm + 0.5 * real(sum_Y2_m ...
            - 2*sum(currentParams.C.means{groupIdx}.*XY{groupIdx}',2));
        for yIdx = 1:yDims(groupIdx)
            phi_bm(yIdx) = phi_bm(yIdx) + 0.5 * real( ...
                           + currentParams.C.moments{groupIdx}{yIdx}(:)' * XX{groupIdx}(:) ...
                           );
        end
        currentParams.phi.b(currGroup(1):currGroup(2)) = phi_bm;
        currentParams.phi.mean(currGroup(1):currGroup(2)) ...
            = currentParams.phi.a ./ phi_bm;
    end
    % Set minimum private variance
    currentParams.phi.mean = min(1./varFloor, currentParams.phi.mean);
    currentParams.phi.b = currentParams.phi.a ./ currentParams.phi.mean;
    
    % ======================
    % Compute lower bound
    % ======================
    
    if ~isnan(lbi)
        lbold = lbi;
    end
    
    if getLB
        % Likelihood term
        log_phi = digammaa_phi - log(currentParams.phi.b); % (yDim x 1) array
        lbi = const_lik + (NT/2) * sum(log_phi) ...
            - sum(currentParams.phi.mean .* (currentParams.phi.b - prior.phi.b));

        % X KL term
        lbi = lbi + (xDim*NT)/2 + res.lb_gp + 0.5*lbterms.logdet_Sx_post;

        log_alpha = nan(numGroups,xDim);
        for groupIdx = 1:numGroups
            log_alpha(groupIdx,:) ...
                = digammaa_alpha(groupIdx) - log(currentParams.alpha.b(groupIdx,:));        
        end

        % C KL term
        lbi = lbi + 0.5*logdetC + 0.5*yDim*xDim;
        for groupIdx = 1:numGroups
            lbi = lbi + (yDims(groupIdx)/2)*sum(log_alpha(groupIdx,:));
            for yIdx = 1:yDims(groupIdx)
                lbi = lbi - 0.5*diag(currentParams.C.moments{groupIdx}{yIdx})'*currentParams.alpha.mean(groupIdx,:)'; 
            end
        end

        % alpha KL term
        lbi = lbi + (numGroups*xDim) * (alogb_alpha - loggammaa_alpha_prior);
        for groupIdx = 1:numGroups
            lbi = lbi + sum( ...
                -currentParams.alpha.a(groupIdx) .* log(currentParams.alpha.b(groupIdx,:)) ...
                - prior.alpha.b .* currentParams.alpha.mean(groupIdx,:) ...
                + (prior.alpha.a - currentParams.alpha.a(groupIdx)) .* (log_alpha(groupIdx,:)) ...
                ) ...
                + xDim * loggammaa_alpha_post(groupIdx) ...
                + xDim * currentParams.alpha.a(groupIdx);
        end

        % phi KL term
        lbi = lbi + yDim * (alogb_phi + loggammaa_phi_post - loggammaa_phi_prior + currentParams.phi.a) ...
            + sum( ...
                  -currentParams.phi.a .* log(currentParams.phi.b) ...
                  - prior.phi.b .* currentParams.phi.mean ...
                  + (prior.phi.a - currentParams.phi.a) .* (digammaa_phi - log(currentParams.phi.b)) ...
              );

        % d KL term
        lbi = lbi + yDim/2 + (yDim/2)*log(prior.d.beta) + 0.5*sum(log(currentParams.d.cov)) ...
            - 0.5 * prior.d.beta * sum(dd);
    else
        lbi = nan;    
    end
    
    % Finish tracking EM iteration time
    tEnd    = toc;
    iterTime = [iterTime tEnd];
    
    % Check stopping conditions or errors
    if verbose
        if getLB
            fprintf('EM iteration %3d of %d        lb %f\r', i, maxIters, lbi);
        else
            fprintf('EM iteration %3d of %d\r', i, maxIters);
        end
    end
    
    lb = [lb lbi];
    
    if i <= 2
        lbbase = lbi;
    elseif (lbi < lbold)
        flags.DecreasingLowerBound = 1;
        disp('Error: Decreasing lower bound');
    elseif ((lbi - lbbase) < (1+tol)*(lbold - lbbase))
        break;
    end
    
end

if (length(lb) < maxIters) && (xDim > 0)
    flags.Convergence = 1;
end

if verbose
    if flags.Convergence == 1
        fprintf('Lower bound converged after %d iterations.\n', length(lb));
    elseif (length(lb) < maxIters) && (xDim <= 0)
        fprintf('Fitting stopped because no significant latent dimensions remain.\n');
    else
        fprintf('Fitting stopped after maxIters (%d) was reached.\n', maxIters);
    end 
end

if any(currentParams.phi.mean == 1./varFloor)
    flags.PrivateVarianceFloor = 1; 
end

% Save tracked parameters
trackedParams.lb = lb;
trackedParams.iterTime = iterTime;
trackedParams.gp_params.Ds = gp_params.Ds;
switch currentParams.covType
    case {'rbf', 'exp'}
        trackedParams.gp_params.gams = gp_params.gams;
    case {'sg', 'expcos'}
        trackedParams.gp_params.gams = gp_params.gams;
        trackedParams.gp_params.nus = gp_params.nus;
end
trackedParams.alphas = alphas;

% Remove estimates of X.cov, if desired
if ~saveXcov
    currentParams = rmfield(currentParams,'X');
end
