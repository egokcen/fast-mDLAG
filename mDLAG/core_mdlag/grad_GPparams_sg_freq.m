function [f,df] = grad_GPparams_sg_freq(p,precomp,const)
%
% [f, df] = grad_GPparams_sg_freq(p, precomp, const)  
%
% Description: Gradient computation for spectral Gaussian GP parameter 
%              optimization, using a frequency domain approximation. This
%              function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma; nu ]
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
%     03 Sep 2023 -- Initial full revision.
%     07 Aug 2024 -- Updated gradient and lower bound with much more
%                    efficient computation of XX.
 
df = zeros(size(p));
f = 0;
gamma = exp(p(1)) + const.minGamma;
nu = p(2);
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;
    
    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';
    freqsmin = freqs-nu;
    freqspls = freqs+nu;
    sqfreqsmin = (2*pi.*(freqsmin)).^2;
    sqfreqspls = (2*pi.*(freqspls)).^2;
    
    % Construct prior spectrum and inverse
    sqexpmin = exp(-0.5*sqfreqsmin./gamma);
    sqexppls = exp(-0.5*sqfreqspls./gamma);
    S = (1-const.eps).*sqrt(pi/(2*gamma)) .* (sqexpmin + sqexppls) + const.eps; % (xDim x 1) vector of diagonal elements
    S_inv = 1./S;
    logdet_S = sum(log(S));
    
    % gamma
    dS_dgamma = (1-const.eps) * sqrt(pi/8).*( ...
        -gamma^(-3/2) .* (sqexpmin + sqexppls) ...
        + gamma^(-5/2) .* (sqfreqsmin.*sqexpmin + sqfreqspls.*sqexppls));   % (xDim x 1) vector of diagonal elements

    Sinv_dSdgamma = S_inv .* dS_dgamma;
    XX_Sinv = precomp.Tu(j).XX .* S_inv;
    df(1) = df(1) - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdgamma);

    % nu
    dS_dnu = (1-const.eps) * sqrt(8*pi^5/gamma^3) ...
        .* (freqsmin.*sqexpmin - freqspls.*sqexppls);   % (xDim x 1) vector of diagonal elements
    Sinv_dSdnu = S_inv .* dS_dnu;
    df(2) = df(2) - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdnu);

    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_S - 0.5 * sum(XX_Sinv);
end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
