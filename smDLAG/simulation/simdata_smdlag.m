function [seq, params] = simdata_smdlag(N, T, binWidth, yDims, xDim, hyperparams, snr, tau, eps, D, Z, varargin)
% 
% [seq, params] = simdata_mdlag(N, T, binWidth, yDims, xDim, hyperparams, snr, tau, eps, D, Z, ...)
%
% Description: Generate simulated data according to a (sparse multi-group) 
%              Delayed Latents Across Groups (smDLAG) model.
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
%     tau       -- (1 x xDim) array; Gaussian process (GP) timescales (ms)
%     eps       -- (1 x xDim) array; GP noise variances
%     D         -- (numGroups x xDim) array; delay matrix (ms)
%     Z         -- (xDim x 1) cell array; Z{j} is a (1 x S_j) array
%                  with S_j inducing point locations for latent j (in ms)
%
%     Optional:
%
%     inducefield -- string; Name of inducing point data field in seq
%                    (default: 'wsm')
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
%                (inducefield) -- (xDim x 1) cell array; element j is 
%                                 a (1 x S_j) array of the inducing 
%                                 points for latent j
%                (latentfield) -- (numGroups*xDim x T) array; 
%                                 delayed latent sequences
%                (obsfield) -- (yDim x T) array; observation sequence
%     params  -- Structure containing smDLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
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
%                    gamma   -- (1 x xDim) array; GP timescales in units of 
%                               time are given by 'binWidth ./ sqrt(gamma)'                                                    
%                    eps     -- (1 x xDim) array; GP noise variances
%                    D       -- (numGroups x xDim) array; delays from 
%                               latents to observed variables. NOTE: Delays
%                               are reported as (real-valued) number of 
%                               time-steps.
%                    Z       -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                               array with S_j inducing point locations for
%                               latent j. Units are in time-steps.
%                    xDim    -- int; number of latent variables
%                    yDims   -- (1 x numGroups) array; dimensionalities of 
%                               each observed group
%
% Author: Evren Gokcen
%
% Revision History:
%     28 Nov 2022 -- Initial full revision.

inducefield = 'wsm';
latentfield = 'xsm';
obsfield = 'y';
assignopts(who, varargin);

% Generate smDLAG model parameters
params = generate_params_smdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                tau, eps, D, Z);

% Generate inducing point sequences                            
seq = generate_W_smdlag(params, N, 'inducefield', inducefield);

% Generate latent sequences
seq = generate_XgivenW_smdlag(seq, params, T, 'latentfield', latentfield, ...
                              'inducefield', inducefield);
                   
% Generate observed sequences
seq = generate_obs_smdlag(seq, params, 'latentfield', latentfield, ...
                          'obsfield', obsfield);        