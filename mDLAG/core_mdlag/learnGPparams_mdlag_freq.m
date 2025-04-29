function res = learnGPparams_mdlag_freq(seq, params, constraints, CPhiC, Y0, varargin)
%
% res = learnGPparams_mdlag_freq(seq, params, constraints, CPhiC, Y0, ...)
%
% Description: Update parameters of GP state model given inferred
%              latent states, using a frequency domain approximation.
%
% Arguments:
%
%     Required:
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
%         X.cov      -- data structure whose jth entry, corresponding
%                       to a group of trials of the same length, has fields
%                         T     -- int; number of time steps for this
%                                  trial group
%                         Vsm   -- (xDim*numGroups x xDim*numGroups x T)
%                                  array; posterior covariance at each 
%                                  timepoint
%                         VsmGP -- (numGroups*T x numGroups*T x xDim) 
%                                  array; posterior covariance of each GP
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
%     constraints -- Structure containing constraints on GP parameters:
%                      maxDelay -- float; maximum delay magnitude, in units
%                                  of time steps
%                      minGamma -- float; minimum gamma value 
%                                  (binWidth^2/tau^2), unitless
%     CPhiC      -- (1 x numGroups) cell array; CPhiC{m} is a (xDim x xDim)
%                   array containing the expectation of C_m'*Phi_m*C_m
%     Y0         -- (1 x length(Tu)) cell array; Y0{j} contains the 
%                   mean-subtracted observations for all trials of
%                   length Tu(j)
%
%     Optional:
%
%     MAXITERS -- int; maximum number of line searches (if >0), 
%                 maximum number of function evaluations (if <0), 
%                 for minimize.m (default:-10)
%     learnDelays -- logical; set true to learn delay parameters;
%                    otherwise, delays will remain fixed at their initial
%                    value (default: true)
%     learnTimescales -- logical; set true to learn timescale parameters;
%                        otherwise, timescales will remain fixed at their 
%                        initial value. (default: true)
%                        NOTE: Only 'rbf' and 'exp' kernels are currently
%                              supported.
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters. 
%            Also includes a term used to compute the lower bound.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Aug 2023 -- Initial full revision. 
%     03 Sep 2023 -- Added spectral Gaussian compatibility.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.
%     18 Jun 2024 -- Added learnTimescales option for rbf, exp kernels.
%     29 Aug 2024 -- Added Y0 as an argument.

MAXITERS = -10; % for minimize.m
learnDelays = true;
learnTimescales = true;
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);

precomp = makePrecomp_GPtimescales_freq(seq,params);

switch params.covType
    case 'rbf'

        % ======================
        % Timescale parameters
        % ======================
        
        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        
        % Value of lower bound component involving GP params (excluding delays)
        res.lb_gp = 0;
        for i = 1:xDim
            const = [];
            
            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;

            init_gam = log(params.gamma(i) - constraints.minGamma);
            
            if learnTimescales
                
                % This does the heavy lifting
                [res_gam, f_gp, ~] = minimize(init_gam, ...
                                              'grad_GPparams_rbf_freq', ...
                                              MAXITERS, ...
                                              precomp(i), ...
                                              const);
                
                gamma(i) = exp(res_gam) + constraints.minGamma;
            else
                [f_gp, ~] ...
                    = grad_GPparams_rbf_freq(init_gam, precomp(i), const);
                gamma(i) =  exp(init_gam) + constraints.minGamma;
            end

            res.lb_gp = res.lb_gp + -f_gp(end);
            
        end
        
        res.gamma = gamma;

    case 'sg'

        % ================================
        % Timescale and center frequency
        % ================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        nu = zeros(1,xDim);

        % Value of lower bound component involving GP params
        res.lb_gp = 0;
        for i = 1:xDim
            const = [];

            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;

            init_gam = log(params.gamma(i) - constraints.minGamma);
            init_nu = params.nu(i);
            init_p = [init_gam; init_nu];

            % This does the heavy lifting
            [res_p, f_gp, ~] = minimize(init_p, ...
                                        'grad_GPparams_sg_freq', ...
                                        MAXITERS, ...
                                        precomp(i), ...
                                        const);
            gamma(i) = exp(res_p(1)) + constraints.minGamma;
            nu(i) = res_p(2);
            res.lb_gp = res.lb_gp + -f_gp(end);

        end

        res.gamma = gamma;
        res.nu = nu;

    case 'exp'

        % ======================
        % Timescale parameters
        % ======================
        
        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        
        % Value of lower bound component involving GP params (excluding delays)
        res.lb_gp = 0; 
        for i = 1:xDim
            const = [];
            
            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;
            
            init_gam = log(params.gamma(i) - constraints.minGamma);

            if learnTimescales
            
                % This does the heavy lifting
                [res_gam, f_gp, ~] = minimize(init_gam, ...
                                              'grad_GPparams_exp_freq', ...
                                              MAXITERS, ...
                                              precomp(i), ...
                                              const);
                
                gamma(i) = exp(res_gam) + constraints.minGamma;
            else
                [f_gp, ~] ...
                    = grad_GPparams_exp_freq(init_gam, precomp(i), const);
                gamma(i) =  exp(init_gam) + constraints.minGamma;
            end
            res.lb_gp = res.lb_gp + -f_gp(end);
            
        end
        
        res.gamma = gamma;

    case 'expcos'

        % ================================
        % Timescale and center frequency
        % ================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        nu = zeros(1,xDim);

        % Value of lower bound component involving GP params
        res.lb_gp = 0;
        for i = 1:xDim
            const = [];

            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;

            init_gam = log(params.gamma(i) - constraints.minGamma);
            init_nu = params.nu(i);
            init_p = [init_gam; init_nu];

            % This does the heavy lifting
            [res_p, f_gp, ~] = minimize(init_p, ...
                                        'grad_GPparams_expcos_freq', ...
                                        MAXITERS, ...
                                        precomp(i), ...
                                        const);
            gamma(i) = exp(res_p(1)) + constraints.minGamma;
            nu(i) = res_p(2);
            res.lb_gp = res.lb_gp + -f_gp(end);

        end

        res.gamma = gamma;
        res.nu = nu;
end

% ======================
% Delay parameters
% ======================

if learnDelays && numGroups > 1

    precomp = makePrecomp_GPdelays_freq(seq,params,CPhiC,Y0);

    D = zeros(numGroups, xDim);

    % We don't include delays to the first group in the optimization
    const = [];
    const.maxDelay = constraints.maxDelay;
    for groupIdx = 2:numGroups

        init_delay = params.D(groupIdx,:).';
        init_delay = log(constraints.maxDelay + init_delay) ...
                   - log(constraints.maxDelay - init_delay);

        % This does the heavy lifting
        [res_delay,~,~] = minimize(init_delay, ...
                                    'grad_GPdelays_freq', ...
                                    MAXITERS, ...
                                    precomp(groupIdx), ...
                                    const);
        D(groupIdx,:) = constraints.maxDelay.*tanh(res_delay./2);

    end
    
    res.D = D;
else
    res.D = params.D;
end
