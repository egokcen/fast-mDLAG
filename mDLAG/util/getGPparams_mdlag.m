function gp_params = getGPparams_mdlag(params, binWidth)
%
% gp_params = getGPparams_mdlag(params, binWidth)
%
% Description: GP parameters can be found in params, but they
%              are easier to interpret when given in units of time. This
%              function gets GP parameters from params and 
%              converts them into the units of time corresponding to 
%              binWidth. binWidth should match the binWidth to which the 
%              model was fitted.
%
% Arguments: 
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the relevant fields
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
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%
% Outputs:
%     
%    gp_params -- structure containing mDLAG GP parameters, converted into
%                 units of time (or 1/time).
%                 For 'rbf' or 'exp':
%                     D   -- (numGroups x xDim) array; delays from latents
%                            to observed variables
%                     tau -- (1 x xDim) array; GP timescales  
%                 For 'sg' or 'expcos':
%                     D   -- (numGroups x xDim) array; delays from latents
%                            to observed variables
%                     tau -- (1 x xDim) array; GP timescales  
%                     nu  -- (1 x xDim) array; center frequencies
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Oct 2022 -- Initial full revision.
%     01 Sep 2023 -- Overhauled for better GP kernel modularity.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

switch params.covType

    case 'rbf'
        % Convert delays from bins to units of time
        gp_params.D = binWidth .* params.D;
        
        % Convert timescales to units of time
        gp_params.tau = binWidth ./ sqrt(params.gamma);

    case 'sg'
        % Convert delays from bins to units of time
        gp_params.D = binWidth .* params.D;
        
        % Convert timescales to units of time
        gp_params.tau = binWidth ./ sqrt(params.gamma);

        % Convert to units of 1/time
        gp_params.nu = params.nu ./ binWidth;

    case 'exp'
        % Convert delays from bins to units of time
        gp_params.D = binWidth .* params.D;
        
        % Convert timescales to units of time
        gp_params.tau = binWidth ./ params.gamma;

    case 'expcos'
        % Convert delays from bins to units of time
        gp_params.D = binWidth .* params.D;
        
        % Convert timescales to units of time
        gp_params.tau = binWidth ./ params.gamma;

        % Convert to units of 1/time
        gp_params.nu = params.nu ./ binWidth;

end
