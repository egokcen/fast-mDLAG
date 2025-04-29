function res = learnGPparams_smdlag(seq, params, constraints, varargin)
%
% res = learnGPparams_smdlag(seq, params, constraints, ...)
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
%                     wsm          -- (xDim x 1) cell array; element j is a
%                                     (1 x S_j) array of the inducing 
%                                     points for latent j
% 
%     params -- Structure containing smDLAG model parameters.
%               Contains the fields
%         covType    -- string; type of GP covariance (e.g., 'rbf')
%         gamma      -- (1 x xDim) array; GP timescales in units of 
%                       time are given by 'binWidth ./ sqrt(gamma)'                                                    
%         eps        -- (1 x xDim) array; GP noise variances
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         Z          -- (xDim x 1) cell array; Z{j} is a (1 x S_j) 
%                       array with S_j inducing point locations for
%                       latent j. Units are in time-steps.
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
%         W          -- data structure whose jth entry, corresponding
%                       to a group of trials of the same length, has fields
%                            T      -- int; number of time steps for this
%                                      trial group
%                            cov    -- (S x S) posterior covariance of the  
%                                      inducing points for all latents.
%                            moment -- (S x S) posterior second moment of 
%                                      the inducing points for all latents.
%     constraints -- Structure containing constraints on GP parameters:
%                      maxDelay -- float; maximum delay magnitude, in units
%                                  of time steps
%                      minGamma -- float; minimum gamma value 
%                                  (binWidth^2/tau^2), unitless
%                      learnInducingLocs -- logical; set true to learn 
%                          inducing point locations; otherwise, inducing 
%                          point locations will remain fixed at their 
%                          initial value
%
%     Optional:
%
%     MAXITERS -- int; maximum number of line searches (if >0), 
%                 maximum number of function evaluations (if <0), 
%                 for minimize.m (default:-10)
%     verbose  -- logical that specifies whether to display status messages
%                 (default: false)
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters:
%            gamma, D, Z. Also includes the number of iterations and
%            value of the cost function after updating these parameters via
%            gradient descent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     05 Dec 2022 -- Initial full revision.
%     11 Aug 2024 -- Overhauled with makePrecomp_GP, grad_GPparams to 
%                    revise several supoptimal computations.

MAXITERS  = -10; % for minimize.m
verbose   = false;
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);

switch params.covType
    case 'rbf'
        % If there's more than one type of parameter, put them in the
        % second row of oldParams.
        oldParams.gamma = params.gamma;
        oldParams.D = params.D;
        oldParams.Z = params.Z;
        fname     = 'grad_GPparams';
    % Optional: Insert other covariance functions here    
end

% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);

% Perform precomputations
precomp = makePrecomp_GP(seq,params);

% Loop once for each state dimension (each GP)
currentParams = params;
for i = 1:xDim
    const = [];
    
    % We're not learning GP noise variance, for now
    const.eps = params.eps(i);
    const.minGamma = constraints.minGamma;
    const.maxDelay = constraints.maxDelay;
    const.learnInducingLocs = constraints.learnInducingLocs;
    
    switch fname                
        case 'grad_GPparams'
            % Change of variables for constrained optimization
            init_gam = log(oldParams.gamma(i) - constraints.minGamma);
            % We don't include delays to the first group in the optimization
            init_delay = reshape(oldParams.D(2:end,i), numGroups-1, 1);
            init_delay = log(constraints.maxDelay + init_delay) ...
                       - log(constraints.maxDelay - init_delay);
            init_Z = oldParams.Z{i}';
            init_p = [init_gam; init_delay; init_Z];
            
            if i > 1
                for j = 1:length(Tu)
                    % Update precomputed GP covariances
                    T = Tu(j);
                    precomp(i).Tu(j).KxwKw_inv ...
                        = precomp(i-1).Tu(j).KxwKw_inv;
                    Kxw = construct_Kxw_smdlag(currentParams, T);
                    Kw = construct_Kw_smdlag(currentParams);
                    try
                        precomp(i).Tu(j).KxwKw_inv{i-1} ...
                            = Kxw{i-1} * invChol_mex(Kw{i-1});
                    catch
                        precomp(i).Tu(j).KxwKw_inv{i-1} ...
                            = Kxw{i-1} / Kw{i-1};
                    end
                    
                    % Update WWKwKwx, another precomputed intermediate term
                    for xIdx = 1:(i-1)
                        precomp(i).Tu(j).WWKwKwxCPhiC{xIdx} ...
                            = precomp(i).Tu(j).Wmoment{xIdx} ...
                            * precomp(i).Tu(j).KxwKw_inv{xIdx}';
                    end
                end
            end
    end   
    
    % This does the heavy lifting
    [res_p, f_gp, grad_iters] =...
        minimize(init_p, fname, MAXITERS, precomp(i), const);
    
    switch params.covType
        case 'rbf'
            switch fname                
                case 'grad_GPparams'
                    currentParams.gamma(i) = exp(res_p(1)) + constraints.minGamma;
                    currentParams.D(2:end,i) = constraints.maxDelay.*tanh(res_p(2:numGroups)./2);
                    currentParams.Z{i} = res_p(numGroups+1:end)';
            end        
    end    
    
    if verbose
        fprintf('\nConverged p; xDim:%d, p:%s', i, mat2str(res_p, 3));
    end
end

res.D = currentParams.D;
res.gamma = currentParams.gamma;
res.Z = currentParams.Z;

end

