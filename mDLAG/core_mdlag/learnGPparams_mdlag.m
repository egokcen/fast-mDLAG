function res = learnGPparams_mdlag(seq, params, constraints, varargin)
%
% res = learnGPparams_mdlag(seq, params, constraints, ...)
%
% Description: Update parameters of GP state model given inferred
%              latent states.
%              NOTE: Learning GP noise variance (eps) is currently
%                    unsupported.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%                     xsm          -- ((numGroups*xDim) x T) array; 
%                                     posterior mean at each timepoint
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
%
%     Optional:
%
%     MAXITERS -- int; maximum number of line searches (if >0), 
%                 maximum number of function evaluations (if <0), 
%                 for minimize.m (default:-10)
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters.
%            Also includes the value of the cost function after updating 
%            these parameters via gradient descent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Oct 2022 -- Initial full revision.
%     02 Sep 2023 -- Overhauled to better to handle alernate GP kernels.
%     08 Mar 2024 -- Added exponential, exponential-cosine kernels.

MAXITERS = -10;   % for minimize.m
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);


precomp = makePrecomp_GP(seq,params);

switch params.covType
    case 'rbf'
        % ====================================
        % Timescale and time delay parameters
        % ====================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        D = zeros(numGroups, xDim);

        % Value of lower bound component involving GP params
        res.lb_gp = 0;

        for i = 1:xDim
            const = [];

            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;
            const.maxDelay = constraints.maxDelay;

            init_gam = log(params.gamma(i) - constraints.minGamma);
            % We don't include delays to the first group in the optimization
            init_delay = reshape(params.D(2:end,i), numGroups-1, 1);
            init_delay = log(constraints.maxDelay + init_delay) ...
                       - log(constraints.maxDelay - init_delay);
            init_p = [init_gam; init_delay];

            % This does the heavy lifting
            [res_p, f_gp, ~] = minimize(init_p, ...
                                        'grad_GPparams_rbf', ...
                                        MAXITERS, ...
                                        precomp(i), ...
                                        const);
            gamma(i) = exp(res_p(1)) + constraints.minGamma;
            D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:end)./2);
            res.lb_gp = res.lb_gp + -f_gp(end);

        end

        res.gamma = gamma;
        res.D = D;

    case 'sg'
        % ========================================================
        % Timescale, center frequency, and time delay parameters
        % ========================================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        D = zeros(numGroups, xDim);
        nu = zeros(1,xDim);

        % Value of lower bound component involving GP params
        res.lb_gp = 0;

        for i = 1:xDim
            const = [];
            
            % We're not learning GP noise variance, for now
            const.eps = params.eps(i);
            const.minGamma = constraints.minGamma;
            const.maxDelay = constraints.maxDelay;

            init_gam = log(params.gamma(i) - constraints.minGamma);
            % We don't include delays to the first group in the optimization
            init_delay = reshape(params.D(2:end,i), numGroups-1, 1);
            init_delay = log(constraints.maxDelay + init_delay) ...
                       - log(constraints.maxDelay - init_delay);
            init_nu = params.nu(i);
            init_p = [init_gam; init_delay; init_nu];

            % This does the heavy lifting
            [res_p, f_gp, ~] = minimize(init_p, ...
                                        'grad_GPparams_sg', ...
                                        MAXITERS, ...
                                        precomp(i), ...
                                        const);
            gamma(i) = exp(res_p(1)) + constraints.minGamma;
            D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:numGroups)./2);
            nu(i) = res_p(end);
            res.lb_gp = res.lb_gp + -f_gp(end);

        end

        res.gamma = gamma;
        res.D = D;    
        res.nu = nu;

    case 'exp'

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        D = params.D;

        % Value of lower bound component involving GP params
        res.lb_gp = 0;

        if numGroups > 1 
            % ====================================
            % Timescale and time delay parameters
            % ====================================
    
            for i = 1:xDim
                const = [];
    
                % We're not learning GP noise variance, for now
                const.eps = params.eps(i);
                const.minGamma = constraints.minGamma;
                const.maxDelay = constraints.maxDelay;
    
                init_gam = log(params.gamma(i) - constraints.minGamma);
                % We don't include delays to the first group in the optimization
                init_delay = reshape(params.D(2:end,i), numGroups-1, 1);
                init_delay = log(constraints.maxDelay + init_delay) ...
                           - log(constraints.maxDelay - init_delay);
                init_p = [init_gam; init_delay];
    
                % This does the heavy lifting
                [res_p, f_gp, ~] = minimize(init_p, ...
                                            'grad_GPparams_exp_approx', ...
                                            MAXITERS, ...
                                            precomp(i), ...
                                            const);
                gamma(i) = exp(res_p(1)) + constraints.minGamma;
                D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:end)./2);
                res.lb_gp = res.lb_gp + -f_gp(end);

            end

        else

            % ====================
            % Timescale parameters
            % ====================
    
            for i = 1:xDim
                const = [];
    
                % We're not learning GP noise variance, for now
                const.eps = params.eps(i);
                const.D = params.D(:,i);
                const.minGamma = constraints.minGamma;
    
                init_gam = log(params.gamma(i) - constraints.minGamma);
                init_p = init_gam;
    
                % This does the heavy lifting
                [res_p, f_gp, ~] = minimize(init_p, ...
                                            'grad_GPtimescales_exp', ...
                                            MAXITERS, ...
                                            precomp(i), ...
                                            const);
                gamma(i) = exp(res_p(1)) + constraints.minGamma;
                res.lb_gp = res.lb_gp + -f_gp(end);
    
            end

        end

        res.gamma = gamma;
        res.D = D;

    case 'expcos'
        % ========================================================
        % Timescale, center frequency, and time delay parameters
        % ========================================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        D = params.D; %zeros(numGroups, xDim);
        nu = zeros(1,xDim);

        % Value of lower bound component involving GP params
        res.lb_gp = 0;

        if numGroups > 1
            for i = 1:xDim
                const = [];
                
                % We're not learning GP noise variance, for now
                const.eps = params.eps(i);
                const.minGamma = constraints.minGamma;
                const.maxDelay = constraints.maxDelay;
    
                init_gam = log(params.gamma(i) - constraints.minGamma);
                % We don't include delays to the first group in the optimization
                init_delay = reshape(params.D(2:end,i), numGroups-1, 1);
                init_delay = log(constraints.maxDelay + init_delay) ...
                           - log(constraints.maxDelay - init_delay);
                init_nu = params.nu(i);
                init_p = [init_gam; init_delay; init_nu];
    
                % This does the heavy lifting
                [res_p, f_gp, ~] = minimize(init_p, ...
                                            'grad_GPparams_expcos_approx', ...
                                            MAXITERS, ...
                                            precomp(i), ...
                                            const);
                gamma(i) = exp(res_p(1)) + constraints.minGamma;
                D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:numGroups)./2);
                nu(i) = res_p(end);
                res.lb_gp = res.lb_gp + -f_gp(end);
    
            end

        else
            % =========================================
            % Timescale and center frequency parameters
            % =========================================
    
            for i = 1:xDim
                const = [];
    
                % We're not learning GP noise variance, for now
                const.eps = params.eps(i);
                const.D = params.D(:,i);
                const.minGamma = constraints.minGamma;
    
                init_gam = log(params.gamma(i) - constraints.minGamma);
                init_nu = params.nu(i);
                init_p = [init_gam; init_nu];
    
                % This does the heavy lifting
                [res_p, f_gp, ~] = minimize(init_p, ...
                                            'grad_GPtimescales_expcos', ...
                                            MAXITERS, ...
                                            precomp(i), ...
                                            const);
                gamma(i) = exp(res_p(1)) + constraints.minGamma;
                nu(i) = res_p(end);
                res.lb_gp = res.lb_gp + -f_gp(end);

            end

        end

        res.gamma = gamma;
        res.D = D;    
        res.nu = nu;  

end
