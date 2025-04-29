% startup.m
%
% Description: Define constants used to characterize model selection across
%              methods, as a function of sequence length.
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

dataFile = 'bias_modelselection.mat';  % Name of dataset files

% Dataset size characteristics
numRuns = 40;                     % Number of independent datasets
numWorkers = numRuns;             % Number of parallel cores to use
Ntotal = 50;                      % Total number of trials
Ntrain = Ntotal;                  % Number of train trials
Ntest = Ntotal - Ntrain;          % Number of test trials
train = 1:Ntrain;                 % Define training set
test = Ntrain+1:Ntotal;           % Define test set
Ttotal = 200;                     % Total number of time points per trial
numPartitions = 5;                % Partition data to sweep trials number
Tlist = floor(logspace(1,log10(Ttotal),numPartitions)); % Increasing T sweep values
binWidth = 20;                    % Sample period of ground truth data
yDims = 24;                       % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);        % Total number of groups
xDim = 4;                         % Latent dimensionality
snrPartitions = 4;
runsPerSNR = floor(numRuns / snrPartitions);
snrList = logspace(-2,1,snrPartitions);
snrRunList = reshape(repmat(snrList, runsPerSNR,1),1,[]); % Signal-to-noise ratios on each run
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
% We'll randomly generate parameters within specified limits
lims.tau = [50 50];    % GP timescale range, in ms
lims.eps = [1e-6 1e-6];  % GP noise range
lims.delay = [-20 20];   % Time delay range, in ms
         
%% General constants

TIMECOLOR = 'k';
SPARSECOLOR = '#2F8E00';
FREQCOLOR = '#D35FBC';
TAPERCOLOR = '#FF8F00';
snrColors = {'#FF5555', '#FF8F00', '#D35FBC', '#2F8E00'};