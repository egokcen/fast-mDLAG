% startup.m
%
% Description: Define constants used to characterize the fitting methods as
%              a function of number of observed groups.
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

dataFile = 'scaling_numgroups';  % Name of dataset files

% Dataset size characteristics
numRuns = 20;                     % Number of independent datasets
numWorkers = numRuns;             % Number of parallel cores to use
Ntotal = 200;                     % Total number of trials
Ntrain = 100;                     % Number of train trials
Ntest = Ntotal - Ntrain;          % Number of test trials
train = 1:Ntrain;                 % Define training set
test = Ntrain+1:Ntotal;           % Define test set
Ttotal = 50;                      % Total number of time points per trial
binWidth = 20;                    % Sample period of ground truth data
numGroupsList = [1 2 3 4 6 8 12 24];  % Sweep over number of groups
numPartitions = length(numGroupsList);
yDim = numGroupsList(end);        % Total number of observed dimensions
xDim = 1;                         % Latent dimensionality
covType = 'rbf';                  % Type of covariance kernel
units = 'ms';

% The parameters below will be used to generate one joint group of
% observations

% ARD hyperparameters
MAG = 100;                    % Controls the variance of loading magnitudes
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;

snr = 0.2;  % Signal-to-noise ratio for each group

% We'll randomly generate parameters within specified limits
lims.tau = [100 100];    % GP timescale range, in ms
lims.eps = [1e-6 1e-6];  % GP noise range
lims.delay = [0 20];     % Time delay range, in ms

%% General constants

TIMECOLOR = 'k';
SPARSECOLOR = '#2F8E00';
FREQCOLOR = '#D35FBC';
TAPERCOLOR = '#FF8F00';
