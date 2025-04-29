function seq = freq2time_mdlag(seq,params,varargin)
%
% seq = freq2time_mdlag(seq,params,varargin)
%
% Description: Convert mDLAG latent variables from the frequency domain
%              to the time domain.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId         -- unique trial identifier
%                    T (1 x 1)       -- number of timesteps
%                    y (yDim x T)    -- neural data
%                    xfft (xDim x T) -- latents in frequency domain
%
%     params  -- Structure containing mDLAG model parameters. 
%                Contains the relevant fields
%         D          -- (numGroups x xDim) array; delays from latents to 
%                       observed variables. NOTE: Delays are reported as 
%                       (real-valued) number of time-steps.
%         xDim       -- int; number of latent variables
%         yDims      -- (1 x numGroups) array; dimensionalities of each 
%                       observed group
%
%     Optional:
%
%     infield  -- string; name of frequency domain latents 
%                 (default: 'xfft')
%     outfield -- string; name of time domain latents 
%                 (default: 'xsm')
%
% Outputs:
%
%     seq     -- same as input structure, but with the additional field
%                    xsm (xDim*numGroups x T) -- latents in time domain
%     
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Aug 2023 -- Initial full revision.

infield = 'xfft';
outfield = 'xsm';
assignopts(who,varargin);

numGroups = length(params.yDims);
N = length(seq);

% Initialize output field
for n = 1:N
    seq(n).(outfield) = []; 
end

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

for j = 1:length(Tu)
    T = Tu(j);
    
    % Handle even and odd sequence lengths
    freqs = (-floor(T/2):floor((T-1)/2))./T;
    
    % Process all trials with length T
    nList = find(Tall == T);
    for groupIdx = 1:numGroups
        Q = exp(-1i*2*pi*params.D(groupIdx,:)'*freqs);
        for n = nList
            seq(n).(outfield) = vertcat(seq(n).(outfield), ...
                real(sqrt(T).*ifft(ifftshift(Q.*seq(n).(infield),2),[],2)));
        end
    end
end
