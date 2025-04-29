function [f,df] = grad_GPparams_sg(p,precomp,const)
%
% [f, df] = grad_GPparams_sg(p, precomp, const)  
%
% Description: Gradient computation for spectral Gaussian GP parameter 
%              optimization. This function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma; D(2,i); ....; D(M,i); nu]'
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
%     02 Sep 2023 -- Initial full revision.

  params = precomp.params;
  numGroups = length(params.yDims);
  df = zeros(size(p));
  f = 0;
  gamma = exp(p(1)) + const.minGamma;
  Betaall = [0; p(2:numGroups)];
  Delayall = const.maxDelay.*tanh(Betaall./2);
  nu = p(end);
  for j = 1:length(precomp.Tu)
      T = precomp.Tu(j).T;
      mT = numGroups*T;
      dK_dBetak = zeros(mT);
      
      Delaydif = repmat(Delayall,T,1); 
      Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT);
      deltaT = (precomp.Tdif(1:mT,1:mT) - Delaydif); 
      deltaTsq = deltaT.^2;
      temp1 = (1-const.eps)*exp(-(gamma/2) * deltaTsq);
      temp2 = temp1 .* cos(2*pi*nu * deltaT);
      temp3 = temp1 .* sin(2*pi*nu * deltaT);
      dtemp = gamma * temp2 .* deltaT + 2*pi*nu * temp3;
      
      K = temp2 + const.eps*eye(mT);
      KinvXmoment = K\precomp.Tu(j).Xmoment;
      df_dK = -0.5*(precomp.Tu(j).numTrials*eye(mT) - KinvXmoment)/K;
      % gamma
      dK_dgamma = -0.5*temp2.*deltaTsq;
      df_dgamma = df_dK(:)' * dK_dgamma(:);
      df(1) = df(1) + df_dgamma;
      % nu
      dK_dnu = -2*pi*deltaT .* temp3;
      df_dnu = df_dK(:)' * dK_dnu(:);
      df(end) = df(end) + df_dnu;
      
      dDelayall_dBetaall = (const.maxDelay/2).*(sech(Betaall./2)).^2;
      for k = 2:length(Betaall)
          idx = k:numGroups:numGroups*T;
          dK_dBetak(:,idx) = dtemp(:,idx)*dDelayall_dBetaall(k);
          dK_dBetak(idx,:) = -dtemp(idx,:)*dDelayall_dBetaall(k);
          dK_dBetak(idx,idx) = 0;
          df_dBetak = df_dK(:)' * dK_dBetak(:);
          df(k) = df(k) + df_dBetak;
          dK_dBetak(:,idx) = 0;
          dK_dBetak(idx,:) = 0;
      end
      
      f = f - 0.5*precomp.Tu(j).numTrials*logdet(K) - 0.5*trace(KinvXmoment); 
  end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
