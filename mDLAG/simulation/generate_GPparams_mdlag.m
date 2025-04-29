function gp_params = generate_GPparams_mdlag(xDim, numGroups, covType, lims)
%
% gp_params = generate_GPparams_mdlag(xDim, numGroups, covType, lims)
%
% Description: Randomly generate mDLAG GP parameters, within specified
%              constraints.
%
% Arguments:
%
%     Required:
%
%     xDim      -- int; Number of latent states
%     numGroups -- int; Number of groups
%     covType   -- string; type of GP covariance. The following GP kernels
%                  are currently supported:
%                      'rbf'    -- Radial basis function, or squared
%                                  exponential
%                      'sg'     -- Spectral Gaussian
%                      'exp'    -- Exponential
%                      'expcos' -- Exponential-cosine
%     lims      -- structure containing limits for the kernel-dependent 
%                  GP parameters:
%                      For 'rbf' or 'exp':
%                          tau   -- (1 x 2) array; lower- and upper-bounds 
%                                   of GP timescales, in units of time
%                          eps   -- (1 x 2) array; lower- and upper-bounds
%                                   of GP noise variances
%                          delay -- (1 x 2) array; lower- and upper-bounds
%                                   of delays, in units of time
%                      For 'sg' or 'expcos':
%                          tau   -- (1 x 2) array; lower- and upper-bounds 
%                                   of GP timescales, in units of time
%                          nu    -- (1 x 2) array; lower- and upper-bounds
%                                   of GP center frequencies, in units
%                                   of 1/time
%                          eps   -- (1 x 2) array; lower- and upper-bounds
%                                   of GP noise variances
%                          delay -- (1 x 2) array; lower- and upper-bounds
%                                   of delays, in units of time
%
% Outputs:
%
%     gp_params -- structure containing kernel-dependent GP parameters:
%                      For 'rbf' or 'exp':
%                          tau  -- (1 x xDim) array; Gaussian process (GP) 
%                                  timescales (units of time)
%                          eps  -- (1 x xDim) array; GP noise variances
%                          D    -- (numGroups x xDim) array; delay matrix 
%                                  (units of time)
%                      For 'sg' or 'expcos':
%                          tau  -- (1 x xDim) array; Gaussian process (GP) 
%                                  timescales (units of time)
%                          nu   -- (1 x xDim) array; GP center frequencies
%                                  (units of 1/time)
%                          eps  -- (1 x xDim) array; GP noise variances
%                          D    -- (numGroups x xDim) array; delay matrix 
%                                  (units of time)
% 
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     23 Oct 2022 -- Initial full revision.
%     01 Sep 2023 -- Overhauled for better GP kernel modularity.
%     07 Mar 2024 -- Added exponential and exponential-cosine kernels.

switch covType

    case {'rbf', 'exp'}
        
        % GP timescales and noise variances
        gp_params.tau = lims.tau(1) + (lims.tau(2)-lims.tau(1)).*rand(1,xDim); 
        % Deal with noise variances on a log scale
        gp_params.eps = exp(log(lims.eps(1)) + (log(lims.eps(2))- log(lims.eps(1))).*(rand(1,xDim)));
        
        % Delays
        delays = cell(1,numGroups);
        for groupIdx = 1:numGroups
            if groupIdx <= 1
                % Set delays to first group to 0 time steps
                delays{groupIdx} = zeros(1,xDim); 
            else
                % All other groups have non-zero delays
                delays{groupIdx} = lims.delay(1) + (lims.delay(2)-lims.delay(1)).*rand(1,xDim);
            end
        end
        % Fill in output structure
        gp_params.D = cat(1, delays{:});

    case {'sg', 'expcos'}

        % GP timescales and noise variances
        gp_params.tau = lims.tau(1) + (lims.tau(2)-lims.tau(1)).*rand(1,xDim); 
        % Center frequencies
        gp_params.nu = lims.nu(1) + (lims.nu(2)-lims.nu(1)).*rand(1,xDim);
        % Deal with noise variances on a log scale
        gp_params.eps = exp(log(lims.eps(1)) + (log(lims.eps(2))- log(lims.eps(1))).*(rand(1,xDim)));
        
        % Delays
        delays = cell(1,numGroups);
        for groupIdx = 1:numGroups
            if groupIdx <= 1
                % Set delays to first group to 0 time steps
                delays{groupIdx} = zeros(1,xDim); 
            else
                % All other groups have non-zero delays
                delays{groupIdx} = lims.delay(1) + (lims.delay(2)-lims.delay(1)).*rand(1,xDim);
            end
        end
        % Fill in output structure
        gp_params.D = cat(1, delays{:});

end