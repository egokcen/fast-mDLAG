function S = make_S_mdlag(params, f)
%
% S = make_S_mdlag(params, f)
%
% Constructs the diagonal elements of the GP spectral density matrix for a
% given frequency.
%
% Arguments:
%
%     params -- mDLAG model parameters
%     f      -- float; frequency, in units of (1/timeSteps)
%
% Outputs:
%
%     S      -- (xDim x 1) array; diagonal elements of the GP spectral 
%               density matrix
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     23 Aug 2023 -- Initial full revision.
%     07 Mar 2024 -- Added exponential and exponential-cosine kernels.
%     12 Aug 2024 -- Only return a vector of diagonal elements, not
%                    a full diagonal matrix.

S = nan(params.xDim,1); % The diagonal elements of S

% Fill in S
for j = 1:params.xDim
    switch params.covType
        case 'rbf'
            S(j) = (1 - params.eps(j)) .* sqrt(2*pi/params.gamma(j)) ...
                   .*exp(-0.5*((2*pi.*f).^2)./params.gamma(j)) ...
                   + params.eps(j);

        case 'sg'
            S(j) = (1 - params.eps(j)) .* sqrt(0.5*pi/params.gamma(j)) ...
                   .*( exp(-0.5*((2*pi.*(f - params.nu(j))).^2)./params.gamma(j)) ...
                    + exp(-0.5*((2*pi.*(f + params.nu(j))).^2)./params.gamma(j)) ) ...
                    + params.eps(j);

        case 'exp'
            S(j) = (1 - params.eps(j)) .* (2*params.gamma(j)) ...
                   ./ (params.gamma(j).^2 + (2*pi.*f).^2) ...
                   + params.eps(j);

        case 'expcos'
            S(j) = (1 - params.eps(j)) .* params.gamma(j) ...
                   .* ( 1 ./ ( (2*pi*(f - params.nu(j))).^2 + params.gamma(j)^2 ) ...
                      + 1 ./ ( (2*pi*(f + params.nu(j))).^2 + params.gamma(j)^2 ) ) ...
                   + params.eps(j);
    end     
end
