function [f,df] = grad_GPparams(p,precomp,const)
%
% [f, df] = grad_GPparams(p, precomp, const)  
%
% Description: Gradient computation for GP timescale, delay, and inducing
%              point location optimization. This function is called by 
%              minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ gamma; D(2,i); ....; D(M,i); Z{i}']'
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
%     05 Dec 2022 -- Initial full revision. 
%     11 Aug 2024 -- Overhauled with learnGPparams, makePrecomp_GP to 
%                    revise several supoptimal computations.
 
numGroups = precomp.numGroups;
xIdx = precomp.xIdx;
xDim = length(precomp.Tu(1).Wmoment);
not_xIdx = [1:xIdx-1 xIdx+1:xDim];  % All but the current latent

df = zeros(size(p));
f = 0;
gamma = exp(p(1)) + const.minGamma;
Betaall = [0; p(2:numGroups)];
Delayall = const.maxDelay.*tanh(Betaall./2);
Z = p(numGroups+1:end);
S = length(Z);
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;
    group_block_idxs = get_block_idxs(repmat(T,1,numGroups));

    % Construct Kw and its inverse
    Zdif = repmat(Z, 1, S) - repmat(Z', S, 1);
    deltaZsq = Zdif.^2;
    temp = (1-const.eps)*exp(-0.5*gamma*deltaZsq); 
    Kw = temp + const.eps*eye(S); % (S x S)
    Kw_inv = inv(Kw);
    CPhiYWKw_inv = precomp.Tu(j).CPhiYW * Kw_inv;  % An intermediate term
    
    % Intermediate terms for derivatives
    dKw_dgamma = -0.5 * temp .* deltaZsq;
    if const.learnInducingLocs
        ztemp_Kw = -gamma * temp .* Zdif;     % Will be used for Z derivative
    end
    
    % Construct Kxw
    Tall = repmat(1:T,numGroups,1)'; % (T x m)
    Dall = repmat(Delayall,1,T)'; % (T x m)
    Tdif = Tall(:) - Dall(:); % (mT x mT), M -> T
    deltaT = Tdif - Z'; % (M -> T) x S
    deltaTsq = deltaT.^2;
    Kxw = (1-const.eps)*exp(-0.5*gamma*deltaTsq); % (M -> T) x S
    KxwKw_inv = Kxw * Kw_inv;
    
    % Intermediate terms for derivatives
    dKwx_dgamma = -0.5 * (Kxw .* deltaTsq)';
    dtemp = -gamma * (Kxw .* deltaT)'; % Will be used for  D & Z derivative
    precomp.Tu(j).WWKwKwx{xIdx} ...
        = precomp.Tu(j).Wmoment{xIdx} * KxwKw_inv';
    
    % Derivative of f w.r.t. gamma
    % Terms w.r.t Kwx and Kxw
    df_dgamma = reshape(CPhiYWKw_inv',[],1)' * dKwx_dgamma(:) ...
        + precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum(KxwKw_inv.*dKwx_dgamma',2);
    KwdKwx_dgamma = Kw_inv * dKwx_dgamma; % An intermediate term
    for i = 1:xDim
        df_dgamma = df_dgamma - precomp.Tu(j).CPhiC_big{i}' * sum(precomp.Tu(j).WWKwKwx{i}' .* KwdKwx_dgamma',2);
    end
    df_dgamma = 2*df_dgamma;
    % Terms w.r.t. Kw
    KwdKw_dgamma = Kw_inv * dKw_dgamma;  % An intermediate term
    df_dgamma = df_dgamma ...
        - 0.5 * reshape(precomp.Tu(j).numTrials * Kw_inv - Kw_inv * precomp.Tu(j).Wmoment{xIdx} * Kw_inv,[],1)' * dKw_dgamma(:) ...
        - reshape(KxwKw_inv,[],1)' * reshape(precomp.Tu(j).CPhiYW * KwdKw_dgamma,[],1) ...
        - 0.5 * precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum((KxwKw_inv * dKw_dgamma) .* KxwKw_inv,2);
    for i = 1:xDim
        df_dgamma = df_dgamma ...
            + precomp.Tu(j).CPhiC_big{i}' * sum(KxwKw_inv .* (precomp.Tu(j).WWKwKwx{i}' * (Kw_inv * dKw_dgamma)),2);
    end

    df(1) = df(1) + df_dgamma;
    
    % Derivative of f w.r.t. delays
    dKwx_dBetak = zeros(size(dtemp));
    dDelayall_dBetaall = (const.maxDelay/2).*(sech(Betaall./2)).^2;
    for k = 2:length(Betaall)
        groupBlock = group_block_idxs{k};
        dKwx_dBetak(:,groupBlock(1):groupBlock(2)) ...
            = -dtemp(:,groupBlock(1):groupBlock(2))*dDelayall_dBetaall(k);
        df_dBetak = reshape(CPhiYWKw_inv',[],1)' * dKwx_dBetak(:) ...
            + precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum(KxwKw_inv.*dKwx_dBetak',2);
        KwdKwx_dBetak = Kw_inv * dKwx_dBetak; % An intermediate term
        for i = 1:xDim
            df_dBetak = df_dBetak - precomp.Tu(j).CPhiC_big{i}' * sum(precomp.Tu(j).WWKwKwx{i}' .* KwdKwx_dBetak',2);
        end
        df(k) = df(k) + 2 * df_dBetak;
        dKwx_dBetak(:,groupBlock(1):groupBlock(2)) = 0;
    end
    
    % Derivative of f w.r.t. Z
    if const.learnInducingLocs
        dKw_dZk = zeros(S);
        dKwx_dZk = zeros(size(dtemp));
        for k = 1:S
            % Terms w.r.t. Kxw and Kwx
            dKwx_dZk(k,:) = -dtemp(k,:);
            df_dZk = reshape(CPhiYWKw_inv',[],1)' * dKwx_dZk(:) ...
                + precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum(KxwKw_inv.*dKwx_dZk',2);
            KwdKwx_dZk = Kw_inv * dKwx_dZk; % An intermediate term
            for i = 1:xDim
                df_dZk = df_dZk - precomp.Tu(j).CPhiC_big{i}' * sum(precomp.Tu(j).WWKwKwx{i}' .* KwdKwx_dZk',2);
            end
            df_dZk = 2 * df_dZk;
            dKwx_dZk(k,:) = 0;

            % Terms w.r.t. Kw
            dKw_dZk(:,k) = -ztemp_Kw(:,k);
            dKw_dZk(k,:) = ztemp_Kw(k,:);
            KwdKw_dZk = Kw_inv * dKw_dZk;  % An intermediate term
            df_dZk = df_dZk ...
                - 0.5 * reshape(precomp.Tu(j).numTrials * Kw_inv - Kw_inv * precomp.Tu(j).Wmoment{xIdx} * Kw_inv,[],1)' * dKw_dZk(:) ...
                - reshape(KxwKw_inv,[],1)' * reshape(precomp.Tu(j).CPhiYW * KwdKw_dZk,[],1) ...
                - 0.5 * precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum((KxwKw_inv * dKw_dZk) .* KxwKw_inv,2);
            for i = 1:xDim
                df_dZk = df_dZk ...
                    + precomp.Tu(j).CPhiC_big{i}' * sum(KxwKw_inv .* (precomp.Tu(j).WWKwKwx{i}' * (Kw_inv * dKw_dZk)),2);
            end
            dKw_dZk(:,k) = 0;
            dKw_dZk(k,:) = 0;

            df(numGroups+k) = df(numGroups+k) + df_dZk;
        end
    end

    % Evaluate objective function
    f = f - 0.5 * precomp.Tu(j).numTrials * logdet(Kw) ...
          - 0.5 * Kw_inv(:)' * precomp.Tu(j).Wmoment{xIdx}(:) ...
          + reshape(KxwKw_inv,[],1)' * precomp.Tu(j).CPhiYW(:) ...
          - 0.5 * precomp.Tu(j).numTrials * T * sum(precomp.Tu(j).CPhiC{xIdx}) ...
          + 0.5 * precomp.Tu(j).numTrials * precomp.Tu(j).CPhiC_big{xIdx}' * sum(KxwKw_inv.*Kxw,2) ...
          - 0.5 * precomp.Tu(j).CPhiC_big{xIdx}' * sum(KxwKw_inv .* precomp.Tu(j).WWKwKwx{xIdx}',2);
    for i = not_xIdx
        f = f - precomp.Tu(j).CPhiC_big{i}' * sum(KxwKw_inv .* precomp.Tu(j).WWKwKwx{i}',2);
    end
end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma - minGamma))
df = -df;
