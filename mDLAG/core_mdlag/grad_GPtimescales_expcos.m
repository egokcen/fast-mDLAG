function [f,df] = grad_GPtimescales_expcos(p,precomp,const)
%
% [f, df] = grad_GPtimescales_expcos(p, precomp, const)  
%
% Description: Gradient computation for exponential-cosine GP timescale 
%              and frequency optimization. This function is called by 
%              minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma; nu]'
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
%     08 Mar 2024 -- Initial full revision.

  params = precomp.params;
  numGroups = length(params.yDims);
  df = zeros(size(p));
  f = 0;
  gamma = exp(p(1)) + const.minGamma;
  nu = p(end);
  for j = 1:length(precomp.Tu)
      T = precomp.Tu(j).T;
      mT = numGroups*T;
      
      Delaydif = repmat(const.D,T,1); 
      Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT);
      deltaT = (precomp.Tdif(1:mT,1:mT) - Delaydif); 
      temp1 = (1-const.eps)*exp(-gamma * abs(deltaT));
      temp2 = temp1 .* cos(2*pi*nu * deltaT);
      temp3 = temp1 .* sin(2*pi*nu * deltaT);
      
      K = temp2 + const.eps*eye(mT);
      KinvXmoment = K\precomp.Tu(j).Xmoment;
      df_dK = -0.5*(precomp.Tu(j).numTrials*eye(mT) - KinvXmoment)/K;
      % gamma
      dK_dgamma = -temp2.*abs(deltaT);
      df_dgamma = df_dK(:)' * dK_dgamma(:);
      df(1) = df(1) + df_dgamma;
      % nu
      dK_dnu = -2*pi*deltaT .* temp3;
      df_dnu = df_dK(:)' * dK_dnu(:);
      df(end) = df(end) + df_dnu;
      
      f = f - 0.5*precomp.Tu(j).numTrials*logdet(K) - 0.5*trace(KinvXmoment); 
  end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
