function params = generate_params_mdlag(yDims, xDim, binWidth, ...
                                        hyperparams, snr, covType, ...
                                        gp_params)
% params = generate_params_mdlag(...)
%
% Description: Randomly generate mDLAG model parameters, within 
%              specified constraints. The params output structure is
%              compatible with other mDLAG codepack functions.
%
% Arguments:
%
%     yDims     -- (1 x numGroups) array; List of dimensionalities of
%                  observed data, [y1Dim, y2Dim, ...]
%     xDim      -- int; Dimensionality of latents, X
%     binWidth  -- float; intended spike count bin width or sample period 
%                  (in units of time). Assume uniform sampling.
%     hyperparams -- structure with the following fields:
%                      beta  -- positive float; precision of mean parameter
%                               generative model (Gaussian)
%                      a_phi -- positive float; 'a' shape parameter of
%                               observation precision (phi) generative 
%                               model (Gamma with mean a/b)
%                      b_phi -- positive float; 'b' scale parameter of
%                               observation precision (phi) generative
%                               model (Gamma with mean a/b)
%                      a_alpha -- (numGroups x xDim) array; 'a' shape 
%                                 parameter of alpha parameter generative
%                                 model (Gamma with mean a/b)
%                      b_alpha -- (numGroups x xDim) array; 'b' scale 
%                                 parameter of alpha parameter generative 
%                                 model (Gamma with mean a/b)
%     snr       -- (1 x numGroups) array; List of signal-to-noise ratios,
%                  defined as trace(CC') / sum(1./phi)
%     covType   -- string; type of GP covariance. The following GP kernels
%                  are currently supported:
%                      'rbf'    -- Radial basis function, or squared
%                                  exponential
%                      'sg'     -- Spectral Gaussian
%                      'exp'    -- Exponential
%                      'expcos' -- Exponential-cosine
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
% Outputs:
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
% 
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     27 Sep 2022 -- Initial full revision.
%     01 Sep 2023 -- Overhauled for better GP kernel modularity.
%     07 Mar 2024 -- Added exponential and exponential-cosine kernels.

numGroups = length(yDims);

% Initialize output structure
params = struct('covType', covType, ...
                'Cs', [], ...
                'alphas', nan(numGroups,xDim), ...
                'phis', [], ...
                'ds', [], ...
                'xDim', xDim, ...
                'yDims', yDims ...
               );
params.Cs = cell(1,numGroups);
params.phis = cell(1,numGroups);
params.ds = cell(1,numGroups);

% Gaussian process parameters
switch covType
    
    case 'rbf'
        params.gamma = (binWidth ./ gp_params.tau).^2;
        params.eps = gp_params.eps;
        params.D = gp_params.D ./ binWidth;        

    case 'sg'
        params.gamma = (binWidth ./ gp_params.tau).^2;
        params.nu = gp_params.nu .* binWidth;
        params.eps = gp_params.eps;
        params.D = gp_params.D ./ binWidth;

    case 'exp'
        params.gamma = binWidth ./ gp_params.tau;
        params.eps = gp_params.eps;
        params.D = gp_params.D ./ binWidth;        

    case 'expcos'
        params.gamma = binWidth ./ gp_params.tau;
        params.nu = gp_params.nu .* binWidth;
        params.eps = gp_params.eps;
        params.D = gp_params.D ./ binWidth;
end

% Generate observation model parameters
for groupIdx = 1:numGroups
    
    for xIdx = 1:xDim
        params.alphas(groupIdx,xIdx) = gamrnd(hyperparams.a_alpha(groupIdx,xIdx), 1./hyperparams.b_alpha(groupIdx,xIdx));
    end
    params.ds{groupIdx} = mvnrnd(zeros(1,yDims(groupIdx)), hyperparams.beta^(-1).*eye(yDims(groupIdx)), 1)';
    params.phis{groupIdx} = gamrnd(hyperparams.a_phi, 1./hyperparams.b_phi, yDims(groupIdx), 1);
    params.Cs{groupIdx} = nan(yDims(groupIdx),xDim);
    for xIdx = 1:xDim
        params.Cs{groupIdx}(:,xIdx) = mvnrnd(zeros(1,yDims(groupIdx)), params.alphas(groupIdx,xIdx)^(-1).*eye(yDims(groupIdx)), 1)';
    end
    
    % Enforce the desired signal-to-noise ratios
    varCC = trace(params.Cs{groupIdx} * params.Cs{groupIdx}');
    varNoise_desired = varCC / snr(groupIdx);
    varNoise_current = sum(params.phis{groupIdx}.^(-1));
    params.phis{groupIdx} = params.phis{groupIdx} .* (varNoise_current / varNoise_desired);
 
end