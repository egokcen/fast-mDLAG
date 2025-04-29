% startup.m
%
% Description: Define constants used throughout the Neuropixels results
%              summaries.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Relevant directories

restoredefaultpath;
addpath ./results
addpath(genpath('../mDLAG'));

%% General constants

% Color schemes
TIMECOLOR = 'k';
SPARSECOLOR = '#2F8E00';
FREQCOLOR = '#D35FBC';

%% For extended demonstrations mDLAG-frequency's scaling

% For scaling with the number of time points per trial
numPartitions = 5;
Tlong = 500;  % Longest pseudo-trial length
Tshort = 50;  % Shortest pseudo-trial length
Tlong_list = ceil(logspace(log10(Tshort), log10(Tlong), numPartitions));

% For scaling with the number of groups
groupNumbers = [2 3 5 7];
