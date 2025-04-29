% check_flags_samplingrate.m
%
% Description: Check the status of model fits across all experiments.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Check status flags

% Time domain
Convergence_time = nan(numPartitions, numRuns);           % Fitting convergence
DecreasingLowerBound_time = nan(numPartitions, numRuns);  % Decreasing lower bound
PrivateVarianceFloor_time = nan(numPartitions, numRuns);  % Private variance floor reached

% Frequency domain
Convergence_freq = nan(numPartitions, numRuns);           % Fitting convergence
DecreasingLowerBound_freq = nan(numPartitions, numRuns);  % Decreasing lower bound
PrivateVarianceFloor_freq = nan(numPartitions, numRuns);  % Private variance floor reached

for runIdx = 1:numRuns
    
    for i = 1:numPartitions
        binWidth = binWidthList(i);
        
        % Load fitted models
        results_time = load( ...
            sprintf('%s/run%03d/period%03d_time.mat', resultDir, runIdx, binWidth) ...
        );
        results_freq = load( ...
            sprintf('%s/run%03d/period%03d_freq.mat', resultDir, runIdx, binWidth) ...
        );
        
        % Flags
        Convergence_time(i,runIdx) = results_time.flags.Convergence;
        DecreasingLowerBound_time(i,runIdx) = results_time.flags.DecreasingLowerBound;
        PrivateVarianceFloor_time(i,runIdx) = results_time.flags.PrivateVarianceFloor;

        Convergence_freq(i,runIdx) = results_freq.flags.Convergence;
        DecreasingLowerBound_freq(i,runIdx) = results_freq.flags.DecreasingLowerBound;
        PrivateVarianceFloor_freq(i,runIdx) = results_freq.flags.PrivateVarianceFloor;
        
    end

end

%% Summarize statuses

% Time domain
fprintf('Convergence, time (numRuns x numPartitions):\n');
disp(Convergence_time')
fprintf('DecreasingLowerBound, time (numRuns x numPartitions):\n');
disp(DecreasingLowerBound_time')
fprintf('PrivateVarianceFloor, time (numRuns x numPartitions):\n');
disp(PrivateVarianceFloor_time')

%% Frequency domain

fprintf('Convergence, freq (numRuns x numPartitions):\n');
disp(Convergence_freq')
fprintf('DecreasingLowerBound, freq (numRuns x numPartitions):\n');
disp(DecreasingLowerBound_freq')
fprintf('PrivateVarianceFloor, freq (numRuns x numPartitions):\n');
disp(PrivateVarianceFloor_freq')
