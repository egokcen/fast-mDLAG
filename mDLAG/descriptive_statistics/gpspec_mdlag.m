function s = gpspec_mdlag(params, binWidth, varargin)
%
% s = gpspec_mdlag(params, binWidth, ...)
%
% Description: Compute and (optionally) plot GP spectra over a range of 
%              frequencies specified by maxfreq.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing mDLAG model parameters.
%                Contains the fields
%         covType -- string; type of GP covariance. The 
%                    following GP kernels are currently 
%                    supported:
%                        'rbf'    -- Radial basis function, or 
%                                    squared exponential
%                        'sg'     -- Spectral Gaussian
%                        'exp'    -- Exponential
%                        'expcos' -- Exponential-cosine
%         The following fields depend on GP covariance type:
%             For 'rbf':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)                                                      
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'sg':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ sqrt(gamma)    
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'exp':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma                                                    
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
%             For 'expcos':
%                 gamma -- (1 x xDim) array; GP timescales; convert to time
%                          via binWidth ./ gamma   
%                 nu    -- (1 x xDim) array; center frequencies; convert to
%                          1/time via nu./binWidth 
%                 eps   -- (1 x xDim) array; GP noise variances
%                 D     -- (numGroups x xDim) array; delays from latents to
%                          observed variables; convert to time via 
%                          D.*binWidth
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
%
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
%     Optional:
%
%     showPlot   -- structure containing a variety of plot options. Specify
%                   as empty ([]) to avoid plotting anything.
%                   amplitude   -- logical; plot amplitude functions
%                                  (default: true)
%                   phasedelay  -- logical; plot phase delay functions
%                                  (default: true)
%                   groupdelay  -- logical; plot group delay functions
%                                  (default: false)
%                   cospec      -- logical; plot co-spectra 
%                                  (default: false)
%                   quadspec    -- logical; plot quadrature spectra
%                                  (default: false)
%                   sqcoherency -- logical; plot squared coherency spectra
%                                  (default: false)
%                   tfgain      -- logical; plot transfer function gains
%                                  (default: false)
%                   intspec     -- logical; plot integrated (cross)-spectra
%                                  (default: false)
%     maxfreq    -- float; maximum frequency (in 1/time, where the units
%                   match that of binWidth) to consider when computing 
%                   cross spectra. (default: 0.5/binWidth, i.e., the
%                   maximum freuency without aliasing)
%     stepres    -- float; resolution of the computed cross spectra, in Hz 
%                   (default: 0.1).
%     pdUnits    -- string; units with which to display phase functions.
%                   'deg' for degrees, 's' for time delay in seconds,
%                   'ms' for time delay in milliseconds (default: 'deg').
%     gdUnits    -- string; units with which to display group delay.
%                   's' for time delay in seconds, 'ms' for time delay in 
%                   milliseconds (default: 'ms').
%     binToSec   -- float; conversion factor from the units of binWidth
%                   to seconds (default: 1/1000, i.e., binWidth is assumed
%                   to be in ms by default)
%     normalize -- logical; set true to compute normalized spectral
%                  density functions. (default: true)
%
% Outputs:
%
%     s -- (numGroups x numGroups) cell array; s{i,j} is a (1 x xDim) cell
%          array, and s{i,j}{r} is the complex-valued spectral density
%          function between group i and group j, given by latent r. 
%          (units: magnitude would be in variance per Hz). sa is an 
%          anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     01 Feb 2023 -- Initial full revision.
%     02 Sep 2023 -- Overhauled with anonymous functions and ability
%                    to compute/plot many versions of the complex-valued
%                    spectra.
%     07 Mar 2024 -- Added exponential, exponential-cosine kernels.

showPlot.amplitude = true;
showPlot.phasedelay = true;
showPlot.groupdelay = false;
showPlot.cospec = false;
showPlot.quadspec = false;
showPlot.sqcoherency = false;
showPlot.tfgain = false;
showPlot.intspec = false;
maxfreq = 0.5/binWidth;
stepres = 0.1;
pdUnits = 'deg';
gdUnits = 'ms';
binToSec = 1/1000;
normalize = true;
assignopts(who,varargin);

binToHz = 1/binToSec; % Convert frequency in binWidth units to Hz
maxfreq = maxfreq * binToHz; % Convert maxfreq to Hz

numGroups = length(params.yDims);
xDim = params.xDim;
freqSteps = (-maxfreq:stepres:maxfreq);  % Hz
freqSteps_half = (0:stepres:maxfreq);    % Hz

% Convert GP params to same units as binWidth
gp_params = getGPparams_mdlag(params,binWidth);

% Unit conversions for phase delay
switch pdUnits
    case 'rad'
        pdfactor = 1;
        pdylbl = 'Phase (rad)';
    case 'deg'
        pdfactor = 180 / pi;
        pdylbl = 'Phase (deg)';
    case 's'
        pdfactor = (2*pi*freqSteps).^(-1);
        pdylbl = 'Phase delay (s)';
    case 'ms'
        pdfactor = (2*pi*freqSteps./binToHz).^(-1);
        pdylbl = 'Phase delay (ms)';
end

% Unit conversions for group delay
switch gdUnits
    case 's'
        gdfactor = (2*pi).^(-1);
        gdylbl = 'Group delay (s)';
    case 'ms'
        gdfactor = (2*pi*binToSec).^(-1);
        gdylbl = 'Group delay (ms)';
end

% Compute spectral densities
s = cell(numGroups,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        s{groupIdx1,groupIdx2} = cell(1,xDim);
        if ~normalize
            CC1 = diag(sum(cat(3,params.C.moments{groupIdx1}{:}),3));
            CC2 = diag(sum(cat(3,params.C.moments{groupIdx2}{:}),3));
        end
        for xIdx = 1:xDim
            D = (gp_params.D(groupIdx2,xIdx) - gp_params.D(groupIdx1,xIdx)).*binToSec; % Convert to sec
            switch params.covType
                case 'rbf'
                    tau = gp_params.tau(xIdx).*binToSec; % Convert to sec
                    s{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f) sqrt(2*pi*tau^2) ...
                        .*exp(-0.5*tau^2.*(2*pi.*f).^2) ...
                        .*exp(1i*2*pi*D.*f);
                case 'sg'
                    tau = gp_params.tau(xIdx).*binToSec; % Convert to sec
                    nu = gp_params.nu(xIdx).*binToHz; % Convert to Hz
                    s{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f)sqrt(0.5*pi*tau^2) ...
                        .*( exp(-0.5*tau^2.*(2*pi.*(f - nu)).^2) ...
                          + exp(-0.5*tau^2.*(2*pi.*(f + nu)).^2)) ...
                        .*exp(1i*2*pi*D.*f);
                case 'exp'
                    gamma = 1./(gp_params.tau(xIdx).*binToSec);
                    s{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f) (2*gamma)./((2*pi*f).^2 + gamma^2) ...
                        .*exp(1i*2*pi*D.*f);
                case 'expcos'
                    gamma = 1/(gp_params.tau(xIdx).*binToSec);
                    nu = gp_params.nu(xIdx).*binToHz; % Convert to Hz
                    s{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f) gamma...
                        .* ( 1 ./ ( (2*pi*(f-nu)).^2 + gamma^2 ) ...
                        + 1 ./ ( (2*pi*(f+nu)).^2 + gamma^2 ) ) ...
                        .*exp(1i*2*pi*D.*f);
            end
            if ~normalize
                s{groupIdx1,groupIdx2}{xIdx} = @(t) sqrt(CC1(xIdx))...
                        .* s{groupIdx1,groupIdx2}{xIdx}(t) .* sqrt(CC2(xIdx));
            end
        end
    end
end

% Plotting
if ~isempty(showPlot)

    if showPlot.amplitude
        % Amplitude functions
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    Amp_curr = abs(s{groupIdx1,groupIdx2}{xIdx}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, Amp_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(Amp_curr);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude (var/Hz)');
                    axis([0 maxfreq 0 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.phasedelay
        % Phase functions
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    Phase_curr = angle(s{groupIdx1,groupIdx2}{xIdx}(freqSteps)).*pdfactor;
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, Phase_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(Phase_curr));
                    if maxval == 0
                       maxval = 1; 
                    end
                    line([0 maxfreq], [0 0], 'color', 'k', 'linestyle', '--');
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(pdylbl);
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.groupdelay
        % Group delay functions
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    groupdelay = @(f) gradient(angle(s{groupIdx1,groupIdx2}{xIdx}(f)),stepres).*gdfactor; % Convert to same units as binWidth
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, groupdelay(freqSteps), 'k-', 'linewidth', 1.5);
                    maxval = max(abs(groupdelay(freqSteps)));
                    if maxval == 0
                       maxval = 1; 
                    end
                    line([0 maxfreq], [0 0], 'color', 'k', 'linestyle', '--');
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(gdylbl);
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.cospec
        % Co-spectra
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    cospec_curr = real(s{groupIdx1,groupIdx2}{xIdx}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, cospec_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(cospec_curr));
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Co-spectrum (var/Hz)');
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.quadspec
        % Across-group quadrature spectra
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    quadspec_curr = imag(s{groupIdx1,groupIdx2}{xIdx}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, quadspec_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(quadspec_curr));
                    if maxval == 0
                       maxval = 1; 
                    end
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Quadrature spectrum (var/Hz)');
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.sqcoherency
        % Coherency spectra
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    sqcoher = @(f) abs(s{groupIdx1,groupIdx2}{xIdx}(f)).^2 ...
                        ./ (s{groupIdx1,groupIdx1}{xIdx}(f) ...
                            .* s{groupIdx2,groupIdx2}{xIdx}(f));
                    sqcoher_curr = sqcoher(freqSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, sqcoher_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Sq. coherency');
                    axis([0 maxfreq 0 1.05]);
                end
            end
        end
    end

    if showPlot.tfgain
        % Gains (amplitude of transfer functions)
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    sqcoher = @(f) abs(s{groupIdx1,groupIdx2}{xIdx}(f)).^2 ...
                        ./ (s{groupIdx1,groupIdx1}{xIdx}(f) ...
                            .* s{groupIdx2,groupIdx2}{xIdx}(f));
                    tfgain = @(f) sqrt((sqcoher(f) ...
                        .* s{groupIdx2,groupIdx2}{xIdx}(f)) ...
                        ./ s{groupIdx1,groupIdx1}{xIdx}(f));
                    tfgain_curr = tfgain(freqSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, tfgain_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(tfgain_curr);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(sprintf('Gain, group %d to %d', groupIdx1, groupIdx2));
                    axis([0 maxfreq 0 1.05*maxval]);
                end
            end
        end
    end

    if showPlot.intspec
        % Integrated spectra
        for xIdx = 1:xDim
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    cospec = @(f) real(s{groupIdx1,groupIdx2}{xIdx}(f));
                    intspec = @(f) cumtrapz(f,2*cospec(f));
                    intspec_curr = intspec(freqSteps_half);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps_half, intspec_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Integrated spectrum (var)');
                    axis([0 maxfreq min(intspec_curr) max([max(intspec_curr) 1.05])]);
                end
            end
        end
    end
end