% startup.m
%
% Description: This script defines constants used throughout the circulant
%              data validations.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

dataDir = './data';        % Access data here
resultDir = './results';   % Access results files here

%% Add directories to matlab path

addpath(dataDir);
addpath(resultDir);
addpath('./util');

%% Synthetic dataset parameters

dataFile = sprintf('data_circulant.mat');    % Name of dataset files

% Dataset size characteristics
N = 100;                            % Total number of trials
T = 100;                            % Number of samples per sequence
binWidth = 20;                      % Sample period of ground truth data
yDims = [10 10 10];                 % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);          % Total number of groups
xDim = 7;                           % Latent dimensionality
snr = 0.2*ones(1,numGroups);        % Signal-to-noise ratios of each group
hyperparams.beta = 1;
hyperparams.a_phi = 1;              
hyperparams.b_phi = 1;
MAG = 100;
hyperparams.a_alpha = MAG.*[1 1   1   Inf 1   Inf Inf;
                            1 1   Inf 1   Inf 1   Inf;
                            1 Inf 1   1   Inf Inf 1  ];
hyperparams.b_alpha = MAG.*ones(numGroups,xDim);
sigDimsTrue = (1./hyperparams.a_alpha) > 0;
units = 'ms';

% Gaussian process (GP) parameters
gp_params.tau = [20 60 40 80 30 150 100];   % GP timescales
gp_params.eps = 1e-5.*ones(1,xDim);         % GP noise variances
gp_params.D = [0   0   0   0  0  0  0;      % Latent delay matrix
               8  -14  0   0  0  0  0;
               22  0  -34 46  0  0  0];
         
%% General constants

groupNames = {'A', 'B', 'C'};
groupColors = {'#5599FF', '#FF5555', '#FF8F00'};
GTCOLOR = '#7F7F7F';
TIMECOLOR = 'k';
FREQCOLOR = '#D35FBC';