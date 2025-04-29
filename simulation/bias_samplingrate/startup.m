% startup.m
%
% Description: Define constants used to characterize the fitting methods as
%              a function of sampling rate.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

dataDir = './data';        % Access data here
resultDir = './results';   % Access results files here
perfDir = './performance_evaluation';  % Performance evaluation scripts

%% Add directories to matlab path

addpath(dataDir);
addpath(resultDir);
addpath(perfDir);

%% Synthetic dataset parameters

dataFile = 'bias_samplingrate.mat';  % Name of dataset files

% Dataset size characteristics
numRuns = 10;                     % Number of independent datasets
numWorkers = numRuns;             % Number of parallel cores to use
Ntotal = 100;                     % Total number of trials
Ntrain = 100;                     % Number of train trials
Ntest = Ntotal - Ntrain;          % Number of test trials
train = 1:Ntrain;                 % Define training set
test = Ntrain+1:Ntotal;           % Define test set
Tlength = 1000;                   % Length of each trial, in ms
binWidthList = [160 120 80 40 20 10]; % List of sample periods to consider
Ttotal = floor(Tlength/min(binWidthList)); % Total number of time points per trial
numPartitions = length(binWidthList); % Number of sample periods
yDims = [12 12];                  % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);        % Total number of groups
xDim = 1;                         % Latent dimensionality
snr = 0.2*ones(1,numGroups);      % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
sigDimsTrue = (1./hyperparams.a_alpha) > 0;
covType = 'rbf';                  % Type of covariance kernel
units = 'ms';

% GP parameters
gp_params.tau = 100;    % GP timescale, in ms
gp_params.D = [0; 10];  % Time delay, in ms
gp_params.eps = 1e-6;   % GP noise
         
%% General constants

GTCOLOR = '#7F7F7F';
TIMECOLOR = 'k';
FREQCOLOR = '#D35FBC';
TAPERCOLOR = '#FF8F00';
