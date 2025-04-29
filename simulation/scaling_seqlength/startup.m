% startup.m
%
% Description: Define constants used to characterize the fitting methods
%              as a function of sequence length.
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

dataFile = 'scaling_seqlength.mat';  % Name of dataset files
longDataFile = 'scaling_long.mat';   % Long-trial datasets

% Dataset size characteristics
numRuns = 20;                     % Number of independent datasets
numWorkers = numRuns;             % Number of parallel cores to use
Ntotal = 200;                     % Total number of trials
Ntrain = 100;                     % Number of train trials
Ntest = Ntotal - Ntrain;          % Number of test trials
train = 1:Ntrain;                 % Define training set
test = Ntrain+1:Ntotal;           % Define test set
Ttotal = 500;                     % Total number of time points per trial
numPartitions = 8;                % Partition data to sweep trials number
Tlist = floor(logspace(1,log10(Ttotal),numPartitions)); % Increasing T sweep values
binWidth = 20;                    % Sample period of ground truth data
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
% For largescale demo of mDLAG-frequency
extraPartitions = 4;
Tlong = 5000;                     
Tlong_list = floor(logspace(log10(Ttotal), log10(Tlong),extraPartitions+1));
Tlong_list = Tlong_list(2:end);

% GP parameters
gp_params.tau = 100;    % GP timescale, in ms
gp_params.D = [0; 10];  % Time delay, in ms
gp_params.eps = 1e-6;   % GP noise
         
%% General constants

TIMECOLOR = 'k';
SPARSECOLOR = '#2F8E00';
FREQCOLOR = '#D35FBC';
TAPERCOLOR = '#FF8F00';