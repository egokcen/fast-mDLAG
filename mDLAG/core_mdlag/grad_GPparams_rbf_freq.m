function [f,df] = grad_GPparams_rbf_freq(p,precomp,const)
%
% [f, df] = grad_GPparams_rbf_freq(p, precomp, const)  
%
% Description: Gradient computation for radial basis (squared exponential)
%              function GP parameter optimization, using a frequency domain
%              approximation. This function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma ]
%     precomp    -- structure containing precomputations
%     const      -- structure containing parameters that stay constant
%                   during this optimization
%
% Outputs:
%
%     f          -- value of portion of lower bound that depends on p
%     df         -- gradient of f at p    
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Aug 2023 -- Initial full revision.
%     07 Aug 2024 -- Updated gradient and lower bound with much more
%                    efficient computation of XX.
 
df = zeros(size(p));
f = 0;
gamma = exp(p(1)) + const.minGamma;
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;
    
    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';
    
    % Construct prior spectrum and inverse
    sqexp = (1 - const.eps).*exp(-0.5*((2*pi.*freqs).^2)./gamma);
    S = sqrt(2*pi/gamma) .* sqexp + const.eps; % (T x 1) vector of diagonal elements
    S_inv = 1./S;
    logdet_S = sum(log(S));
    
    dS_dgamma = sqrt(pi/2) ...
        .* ((2*pi.*freqs).^2 .* gamma^(-5/2) - gamma^(-3/2)) ...
        .* sqexp;   % (xDim x 1) vector of diagonal elements
    Sinv_dSdgamma = S_inv .* dS_dgamma;
    XX_Sinv = precomp.Tu(j).XX .* S_inv;
    df(1) = df(1) - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdgamma);
    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_S - 0.5 * sum(XX_Sinv);
end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
