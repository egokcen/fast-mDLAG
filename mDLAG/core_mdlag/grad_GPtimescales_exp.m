function [f,df] = grad_GPtimescales_exp(p,precomp,const)
%
% [f, df] = grad_GPtimescales_exp(p, precomp, const)  
%
% Description: Gradient computation for exponential covariance timescale
%              optimization. This function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma ]'
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
  for j = 1:length(precomp.Tu)
      T = precomp.Tu(j).T;
      mT = numGroups*T;
      
      Delaydif = repmat(const.D,T,1); 
      Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT);
      deltaT = (precomp.Tdif(1:mT,1:mT) - Delaydif); 
      temp = (1-const.eps)*exp(-gamma * abs(deltaT));
      
      K = temp + const.eps*eye(mT);
      KinvXmoment = K\precomp.Tu(j).Xmoment;
      df_dK = -0.5*(precomp.Tu(j).numTrials*eye(mT) - KinvXmoment)/K;
      dK_dgamma = -temp.*abs(deltaT);
      df_dgamma = df_dK(:)' * dK_dgamma(:);
      df(1) = df(1) + df_dgamma;
      
      f = f - 0.5*precomp.Tu(j).numTrials*logdet(K) - 0.5*trace(KinvXmoment); 
  end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
