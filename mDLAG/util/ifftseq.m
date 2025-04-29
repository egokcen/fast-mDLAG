function seq = ifftseq(seq,infield,outfield)
%
% fftseq
%
% Description: Compute the unitary inverse FFT of each sequence in 
%              seq.(infield). Store in seq.(outfield).
%
% Arguments:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId              -- unique trial identifier
%                     T (1 x 1)            -- number of timesteps
%                     (infield) (dim x T)  -- data sequence
%
% Outputs :
%
%     seq      -- same as the input structure, with additional field
%                     (outfield) (dim x T)  -- unitary inverse FFT of 
%                     seq.(infield)
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     04 Sep 2023 -- Initial full revision.

for n = 1:length(seq)
    T = seq(n).T;
    seq(n).(outfield) = sqrt(T).*ifft(ifftshift(seq(n).(infield),2),[],2);
end
