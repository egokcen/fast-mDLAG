function [seq, params] = simdata_mdlag(N, T, binWidth, yDims, xDim, ...
    hyperparams, snr, covType, gp_params, varargin)
% 
% [seq, params] = simdata_mdlag(N, T, binWidth, yDims, xDim, ...
%                               hyperparams, snr, covType, gp_params, ...)
%
% Description: Generate simulated data according to a (multi-group) Delayed
%              Latents Across Groups (mDLAG) model.
%
% Arguments:
%
%     Required:
%
%     N         -- int; number of sequences
%     T         -- int; number of samples per sequence
%     binWidth  -- float; intended spike count bin width or sample period 
%                  (in units of time). Assume uniform sampling.
%     yDims     -- (1 x M) array; List of dimensionalities of
%                  observed data, [y1Dim, y2Dim, ...]
%     xDim      -- int; Dimensionality of latents, X
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
%     Optional:
%
%     latentfield -- string; Name of latent data field in seq 
%                    (default: 'xsm')
%     obsfield    -- string; Name of observation data field in seq 
%                    (default: 'y')
%
% Outputs:
%
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- int; unique trial (sequence) identifier  
%                T -- int; number of timesteps
%                (latentfield) -- (numGroups*xDim x T) array; 
%                                 delayed latent sequences
%                (obsfield) -- (yDim x T) array; observation sequence
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
%                                   'expcos' -- Exponential cosine
%                    The following fields depend on GP covariance type:
%                        For 'rbf' or 'exp':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma)                                                    
%                            eps   -- (1 x xDim) array; GP noise variances
%                            D     -- (numGroups x xDim) array; delays from 
%                                     latents to observed variables; 
%                                     convert to time via D.*binWidth
%                        For 'sg' or 'expcos':
%                            gamma -- (1 x xDim) array; GP timescales;
%                                     convert to time via 
%                                     binWidth ./ sqrt(gamma) 
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
% Author: Evren Gokcen
%
% Revision History:
%     27 Sep 2022 -- Initial full revision.
%     01 Sep 2023 -- Overhauled for better GP kernel modularity.
%     07 Mar 2024 -- Updated documentation to reflect new GP kernels.

latentfield = 'xsm';
obsfield = 'y';
assignopts(who, varargin);

% Generate mDLAG model parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

% Generate latent sequences
seq = generate_latents_mdlag(params, T, N, 'latentfield', latentfield);
                   
% Generate observed sequences
seq = generate_obs_mdlag(seq, params, 'latentfield', latentfield, ...
                         'obsfield', obsfield);        