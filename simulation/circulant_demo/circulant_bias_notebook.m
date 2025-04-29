% circulant_bias_notebook.m
%
% Description: Investigate the effects of the circulant approximation
%              induced by the frequency domain approach.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_circulant;
addpath('./util');

%% ========================================================================
%  1) Visually explore the differences between time-domain and
%     frequency-domain data generation.
% =========================================================================

% Set mDLAG model parameters. We only need to generate one trial for
% demonstration.
randomSeed = 12;              % Seed for the random number generator
N = 1;                        % Total number of trials
T = 25;                       % Number of samples per sequence
binWidth = 20;                % Sample period of ground truth data
yDims = [1];                  % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);    % Total number of groups
xDim = 1;                     % Latent dimensionality
snr = ones(1,numGroups);      % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
gp_params.tau = 100;                 % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params.D = 0;                     % GP time delays

%% ========================================================================
%  1a) Visualize the GP covariance matrix and its circulant (frequency 
%      domain) counterpart
% =========================================================================

% Generate mDLAG parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

% Time domain covariance matrix
K_time = construct_Kx(params, T); K_time = K_time{1};

% Frequency domain covariance matrix
freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
% Construct the spectral density matrix
S_freq = nan(T,1);
for f = 1:T
    S_freq(f) = make_S_mdlag(params,freqs(f));
end
S_freq = diag(S_freq);
% Construct the time domain covariance matrix from the spectral density
K_freq = real(T.*ifft(ifft(ifftshift(ifftshift(S_freq,1),2),[],1)',[],1));

% Plot both covariance matrices
figure;
% Time domain K
subplot(1,2,1);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(K_time));
axis square;
axis off;
title('Time domain-generated K');

% Frequency domain-generated K
subplot(1,2,2);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(K_freq));
axis square;
axis off;
title('Frequency domain-generated K');

%% ========================================================================
% 1b) Visualize the GP spectral density matrix and its approximate diagonal 
%     (frequency domain) counterpart
% =========================================================================

% Generate mDLAG parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

% Time domain covariance matrix
K_time = construct_Kx(params, T); K_time = K_time{1};
S_time = real((1/T).*fft(fft(K_time,[],1)',[],1));

% Frequency domain covariance matrix
freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
% Construct the spectral density matrix
S_freq = nan(T,1);
for f = 1:T
    S_freq(f) = make_S_mdlag(params,freqs(f));
end
S_freq = diag(ifftshift(S_freq));

% Plot both spectral density matrices (on a log scale)
figure;
subplot(1,2,1);
hold on;
colormap('gray');
colorbar
clim([-5 1.1]);
imagesc(flipud(log10(abs(S_time))));
axis square;
axis off;
title('Time domain-generated S');

subplot(1,2,2);
hold on;
colormap('gray');
colorbar
clim([-5 1.1]);
imagesc(flipud(log10(S_freq)));
axis square;
axis off;
title('Frequency domain-generated S');

%% ========================================================================
%  1c) Compare time domain- and frequency domain-generated signals
% =========================================================================

rng(randomSeed);  % Seed the random number generator

% Time domain generation
[seq_time, params] = simdata_mdlag(N, T, binWidth, yDims, xDim, ...
                                   hyperparams, snr, covType, gp_params);
% Plot shifted copies of the latent signal, to demonstrate the boundary
% conditions (or lack thereof)
figure;
hold on;
plot(binWidth.*(-seq_time(1).T:-1), seq_time(1).xsm, ...
     'color', 'k', ...
     'linestyle', ':', ...
     'linewidth', 1.5);
plot(binWidth.*(0:seq_time(1).T-1), seq_time(1).xsm, ...
     'color', 'k', ...
     'linestyle', '-', ...
     'linewidth', 1.5);
plot(binWidth.*(seq_time(1).T:2*seq_time(1).T-1), seq_time(1).xsm, ...
     'color', 'k', ...
     'linestyle', ':', ...
     'linewidth', 1.5);
ylim([-1.5 1.5]);
xlabel('Time (ms)');
ylabel('x');
title('Time domain generation');

% Frequency domain generation
[seq_freq, params] = simdata_mdlag_freq(N, T, binWidth, yDims, xDim, ...
                                        hyperparams, snr, covType, gp_params);
% Plot shifted copies of the latent signal, to demonstrate the boundary
% conditions (or lack thereof)
figure;
hold on;
plot(binWidth.*(-seq_freq(1).T:-1), seq_freq(1).xsm, ...
     'color', 'k', ...
     'linestyle', ':', ...
     'linewidth', 1.5);
plot(binWidth.*(0:seq_freq(1).T-1), seq_freq(1).xsm, ...
     'color', 'k', ...
     'linestyle', '-', ...
     'linewidth', 1.5);
plot(binWidth.*(seq_freq(1).T:2*seq_freq(1).T-1), seq_freq(1).xsm, ...
     'color', 'k', ...
     'linestyle', ':', ...
     'linewidth', 1.5);
ylim([-1.1 1.1]);
xlabel('Time (ms)');
ylabel('x');
title('Frequency domain generation');

%% ========================================================================
%  2) Visually explore the differences between time-domain and
%     frequency-domain inference.
% =========================================================================

% Set mDLAG model parameters. We only need to generate one trial for
% demonstration.
randomSeed = 20;               % Seed for the random number generator
N = 1;                        % Total number of trials
T = 50;                       % Number of samples per sequence
binWidth = 20;                % Sample period of ground truth data
yDims = [24];                 % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);    % Total number of groups
xDim = 1;                     % Latent dimensionality
snr = 10.*ones(1,numGroups);  % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
gp_params.tau = 100;                 % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params.D = 0;                     % GP time delays

%% ========================================================================
% 2a) Compare time domain- and frequency domain-inferred signals, using
%     the ground truth parameters in both cases.
% =========================================================================

rng(randomSeed);  % Seed the random number generator

% Time domain generation
[seq_true, params] = simdata_mdlag(N, T, binWidth, yDims, xDim, ...
                                   hyperparams, snr, covType, gp_params);
% Compute the FFT of the observed data
seq_true = fftseq(seq_true, 'y', 'yfft');

% Modify the structure of params for compatibility with inference functions
params.C.means = params.Cs;
% Fill in the second moment of C
for groupIdx = 1:numGroups
    for yIdx = 1:yDims(groupIdx)
        % Second moment
        params.C.moments{groupIdx}{yIdx} ...
            = zeros(xDim) ...
            + params.C.means{groupIdx}(yIdx,:)' * params.C.means{groupIdx}(yIdx,:);
    end
end
params.phi.mean = params.phis{1};
params.d.mean = params.ds{1};

% Time domain inference
[seq_time,~,~] = inferX(seq_true, params);

% Frequency domain inference
[seq_freq,~,~] = inferX_freq(seq_true, params);
% Convert estimated latents back to the time domain
seq_freq = freq2time_mdlag(seq_freq, params, 'infield', 'xfft', 'outfield', 'xsm');

% Overlay the signals inferred by each approach on top of the ground truth
figure;
hold on;
h1 = plot(binWidth.*(0:seq_true(1).T-1), seq_true(1).xsm, ...
          'color', GTCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
h2 = plot(binWidth.*(0:seq_time(1).T-1), seq_time(1).xsm, ...
          'color', TIMECOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
h3 = plot(binWidth.*(0:seq_freq(1).T-1), seq_freq(1).xsm, ...
          'color', FREQCOLOR, ...
          'linestyle', '-', ...
          'linewidth', 1.5);
ylim([-1.5 1.5]);
xlabel('Time (ms)');
ylabel('x');
legend([h1 h2 h3], {'Ground truth', 'Time domain inference', ...
    'Frequency domain inference'}, 'location', 'southoutside');

%% ========================================================================
% 3) Quantify the error between the GP covariance matrix and its circulant
%    (frequency domain) approximation.
% =========================================================================

% Visualize example covariance matrices for the multi-group case

% Set mDLAG model parameters
T = 25;                       % Number of samples per sequence
binWidth = 20;                % Sample period of ground truth data
yDims = [1 1];                % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);    % Total number of groups
xDim = 1;                     % Latent dimensionality
snr = 10.*ones(1,numGroups);  % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
gp_params.tau = 100;                 % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params.D = [0; 80];               % GP time delays

% Generate mDLAG parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

% Time domain covariance matrix
K_time = construct_Kx(params, T); 
K_time = cell2mat(K_time);

% Frequency domain covariance matrix
freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
% Construct the spectral density matrix
K_freq = cell(numGroups,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        S_freq = nan(T,1);
        for f = 1:T
            S_freq(f) = make_S_mdlag(params,freqs(f)) .* exp(-1i*2*pi*freqs(f).*(params.D(groupIdx2,:) - params.D(groupIdx1,:)) );
        end
        S_freq = diag(S_freq);
        K_freq{groupIdx1,groupIdx2} = real(T.*ifft(ifft(ifftshift(ifftshift(S_freq,1),2),[],1)',[],1));
    end
end
K_freq = cell2mat(K_freq);

% Plot both covariance matrices
figure;
% Time domain K
subplot(1,2,1);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(K_time));
axis square;
axis off;
title('Time domain-generated K');

% Frequency domain-generated K
subplot(1,2,2);
hold on;
colormap('gray');
colorbar
clim([0 1]);
imagesc(flipud(K_freq));
axis square;
axis off;
title('Frequency domain-generated K');

%% ========================================================================
% 3a) Visualize an example phase-shift matrix
% =========================================================================

% Set mDLAG model parameters
T = 25;                       % Number of samples per sequence
binWidth = 20;                % Sample period of ground truth data
yDims = [1 1];                % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);    % Total number of groups
xDim = 1;                     % Latent dimensionality
snr = 10.*ones(1,numGroups);  % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
gp_params.tau = 100;                 % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params.D = [0; 10];               % GP time delays

% Generate mDLAG parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

% Handle even and odd sequence lengths
freqs = (-floor(T/2):floor((T-1)/2))./T;  
% Construct the phase shift matrix
Ph_freq = cell(numGroups,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        Q_freq = nan(T,1);  % CSD phase
        for f = 1:T
            Q_freq(f) = exp(-1i*2*pi*freqs(f).*(params.D(groupIdx2,:) - params.D(groupIdx1,:)) );
        end
        % Phase-shift matrix, in degrees
        Ph_freq{groupIdx1,groupIdx2} = angle(diag(ifftshift(Q_freq))) .* (180/pi);
    end
end

% Plot the phase-shift matrix to one group
figure;
hold on;
colormap('gray');
colorbar
clim([-90 90]);
imagesc(flipud(Ph_freq{1,2}));
axis square;
axis off;
title('Frequency domain phase shift');

%% ========================================================================
% 3b) Approximation error as a function of trial length
% =========================================================================

% Set mDLAG model parameters
Tlist = ceil(logspace(0,3,100));  % Number of samples per sequence
binWidth = 20;                   % Sample period of ground truth data
yDims = [1 1];                   % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);       % Total number of groups
xDim = 1;                        % Latent dimensionality
snr = 10.*ones(1,numGroups);     % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
gp_params.tau = 100;                 % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params.D = [0; 10];               % GP time delays

% Generate mDLAG parameters
params = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                               covType, gp_params);

K_freq_err = nan(1,length(Tlist));
for t = 1:length(Tlist)
    T = Tlist(t);

    % Time domain covariance matrix
    K_time = construct_Kx(params, T); 
    K_time = cell2mat(K_time);
    
    % Frequency domain covariance matrix
    freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
    % Construct the spectral density matrix
    K_freq = cell(numGroups,numGroups);
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            S_freq = nan(T,1);
            for f = 1:T
                S_freq(f) = make_S_mdlag(params,freqs(f)).*exp(-1i*2*pi*freqs(f).*(params.D(groupIdx2,:) - params.D(groupIdx1,:)) );
            end
            S_freq = diag(S_freq);
            K_freq{groupIdx1,groupIdx2} = real(T.*ifft(ifft(ifftshift(ifftshift(S_freq,1),2),[],1)',[],1));
        end
    end
    K_freq = cell2mat(K_freq);
    
    % Compute approximation error
    K_freq_err(t) = norm(K_time - K_freq, 'fro') / norm(K_time, 'fro');
end

% Visualize asymptotic error
figure;
hold on;
plot(Tlist, K_freq_err, ...
     'color', 'k', ...
     'linestyle', '-', ...
     'linewidth', 1.5);
xlabel('Time points per trial');
ylabel('Approximation error');
xscale('log');
ylim([0 12]);
axis square;

%% ========================================================================
% 3c) Approximation error as a function of estimated GP timescale
% =========================================================================

% Set mDLAG model parameters
T = 25;                          % Number of samples per sequence
binWidth = 20;                   % Sample period of ground truth data
yDims = [1 1];                   % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);       % Total number of groups
xDim = 1;                        % Latent dimensionality
snr = 10.*ones(1,numGroups);     % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
% Time domain parameters
gp_params_time.tau = 100;                 % GP timescales
gp_params_time.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params_time.D = [0; 10];               % GP time delays
% Frequency domain parameters
taulist = linspace(0,200,100);            % GP timescales
gp_params_freq.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params_freq.D = [0; 10];               % GP time delays

K_freq_err = nan(1,length(taulist));
for t = 1:length(taulist)
    % Generate time domain parameters
    params_time = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                        covType, gp_params_time);
    % Generate frequency domain parameters
    gp_params_freq.tau = taulist(t);
    params_freq = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                        covType, gp_params_freq);
    
    % Time domain covariance matrix
    K_time = construct_Kx(params_time, T); 
    K_time = cell2mat(K_time);
    
    % Frequency domain covariance matrix
    freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
    % Construct the spectral density matrix
    K_freq = cell(numGroups,numGroups);
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            S_freq = nan(T,1);
            for f = 1:T
                S_freq(f) = make_S_mdlag(params_freq,freqs(f)).*exp(-1i*2*pi*freqs(f).*(params_freq.D(groupIdx2,:) - params_freq.D(groupIdx1,:)) );
            end
            S_freq = diag(S_freq);
            K_freq{groupIdx1,groupIdx2} = real(T.*ifft(ifft(ifftshift(ifftshift(S_freq,1),2),[],1)',[],1));
        end
    end
    K_freq = cell2mat(K_freq);
    
    % Compute approximation error
    K_freq_err(t) = norm(K_time - K_freq, 'fro') / norm(K_time, 'fro');
end

% Visualize asymptotic error
figure;
hold on;
line([gp_params_time.tau gp_params_time.tau], [0.2 1.2], ...
     'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 1.5);
plot(taulist, K_freq_err, ...
     'color', 'k', ...
     'linestyle', '-', ...
     'linewidth', 1.5);
[min_err, min_idx] = min(K_freq_err);
plot(taulist(min_idx), min_err, ...
     'color', '#FF5555', ...
     'linestyle', 'none', ...
     'marker', '.', ...
     'markersize', 10);
xlabel('Estimated GP timescale (ms)');
ylabel('Approximation error');
ylim([0.2 1.2]);
axis square;

%% ========================================================================
% 3d) Approximation error as a function of estimated GP time delay
% =========================================================================

% Set mDLAG model parameters
T = 25;                          % Number of samples per sequence
binWidth = 20;                   % Sample period of ground truth data
yDims = [1 1];                   % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);       % Total number of groups
xDim = 1;                        % Latent dimensionality
snr = 10.*ones(1,numGroups);     % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
covType = 'rbf';  % Type of covariance kernel

% Gaussian process (GP) parameters
% Time domain parameters
gp_params_time.tau = 100;                 % GP timescales
gp_params_time.eps = 1e-5.*ones(1,xDim);  % GP noise variances
gp_params_time.D = [0; 10];               % GP time delays
% Frequency domain parameters
gp_params_freq.tau = 100;                 % GP timescales
gp_params_freq.eps = 1e-5.*ones(1,xDim);  % GP noise variances
Dlist = linspace(0,20,100);               % GP time delays

K_freq_err = nan(1,length(Dlist));
for t = 1:length(Dlist)
    % Generate time domain parameters
    params_time = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                        covType, gp_params_time);
    % Generate frequency domain parameters
    gp_params_freq.D = [0; Dlist(t)];
    params_freq = generate_params_mdlag(yDims, xDim, binWidth, hyperparams, snr, ...
                                        covType, gp_params_freq);
    
    % Time domain covariance matrix
    K_time = construct_Kx(params_time, T); 
    K_time = cell2mat(K_time);
    
    % Frequency domain covariance matrix
    freqs = (-floor(T/2):floor((T-1)/2))./T;  % Handle even and odd sequence lengths
    % Construct the spectral density matrix
    K_freq = cell(numGroups,numGroups);
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            S_freq = nan(T,1);
            for f = 1:T
                S_freq(f) = make_S_mdlag(params_freq,freqs(f)).*exp(-1i*2*pi*freqs(f).*(params_freq.D(groupIdx2,:) - params_freq.D(groupIdx1,:)) );
            end
            S_freq = diag(S_freq);
            K_freq{groupIdx1,groupIdx2} = real(T.*ifft(ifft(ifftshift(ifftshift(S_freq,1),2),[],1)',[],1));
        end
    end
    K_freq = cell2mat(K_freq);
    
    % Compute approximation error
    K_freq_err(t) = norm(K_time - K_freq, 'fro') / norm(K_time, 'fro');
end

% Visualize asymptotic error
figure;
hold on;
line([gp_params_time.D(end) gp_params_time.D(end)], [0.355 0.362], ...
     'color', [0.5 0.5 0.5], 'linestyle', '--', 'linewidth', 1.5);
plot(Dlist, K_freq_err, ...
     'color', 'k', ...
     'linestyle', '-', ...
     'linewidth', 1.5);
[min_err, min_idx] = min(K_freq_err);
plot(Dlist(min_idx), min_err, ...
     'color', '#FF5555', ...
     'linestyle', 'none', ...
     'marker', '.', ...
     'markersize', 10);
xlabel('Estimated GP time delay (ms)');
ylabel('Approximation error');
ylim([0.355 0.362]);
axis square;