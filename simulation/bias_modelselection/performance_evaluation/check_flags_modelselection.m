% check_flags_modelselection.m
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
xDimsRemoved_time = nan(numPartitions, numRuns);          % Number of latents removed

% Inducing points
Convergence_sparse = nan(numPartitions, numRuns);           % Fitting convergence
DecreasingLowerBound_sparse = nan(numPartitions, numRuns);  % Decreasing lower bound
PrivateVarianceFloor_sparse = nan(numPartitions, numRuns);  % Private variance floor reached
xDimsRemoved_sparse = nan(numPartitions, numRuns);          % Number of latents removed

% Frequency domain
Convergence_freq = nan(numPartitions, numRuns);           % Fitting convergence
DecreasingLowerBound_freq = nan(numPartitions, numRuns);  % Decreasing lower bound
PrivateVarianceFloor_freq = nan(numPartitions, numRuns);  % Private variance floor reached
xDimsRemoved_freq = nan(numPartitions, numRuns);          % Number of latents removed

for runIdx = 1:numRuns
    
    for i = 1:numPartitions
        T = Tlist(i);
        
        % Load fitted models
        results_time = load( ...
            sprintf('%s/run%03d/T%04d_time.mat', resultDir, runIdx, T) ...
        );
        results_sparse = load( ...
            sprintf('%s/run%03d/T%04d_sparse.mat', resultDir, runIdx, T) ...
        );
        results_freq = load( ...
            sprintf('%s/run%03d/T%04d_freq.mat', resultDir, runIdx, T) ...
        );
        
        % Flags
        Convergence_time(i,runIdx) = results_time.flags.Convergence;
        DecreasingLowerBound_time(i,runIdx) = results_time.flags.DecreasingLowerBound;
        PrivateVarianceFloor_time(i,runIdx) = results_time.flags.PrivateVarianceFloor;
        xDimsRemoved_time(i,runIdx) = results_time.flags.xDimsRemoved;

        Convergence_sparse(i,runIdx) = results_sparse.flags.Convergence;
        DecreasingLowerBound_sparse(i,runIdx) = results_sparse.flags.DecreasingLowerBound;
        PrivateVarianceFloor_sparse(i,runIdx) = results_sparse.flags.PrivateVarianceFloor;
        xDimsRemoved_sparse(i,runIdx) = results_sparse.flags.xDimsRemoved;

        Convergence_freq(i,runIdx) = results_freq.flags.Convergence;
        DecreasingLowerBound_freq(i,runIdx) = results_freq.flags.DecreasingLowerBound;
        PrivateVarianceFloor_freq(i,runIdx) = results_freq.flags.PrivateVarianceFloor;
        xDimsRemoved_freq(i,runIdx) = results_freq.flags.xDimsRemoved;
        
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
fprintf('xDimsRemoved, time (numRuns x numPartitions):\n');
disp(xDimsRemoved_time')

%% Inducing points

fprintf('Convergence, sparse (numRuns x numPartitions):\n');
disp(Convergence_sparse')
fprintf('DecreasingLowerBound, sparse (numRuns x numPartitions):\n');
disp(DecreasingLowerBound_sparse')
fprintf('PrivateVarianceFloor, sparse (numRuns x numPartitions):\n');
disp(PrivateVarianceFloor_sparse')
fprintf('xDimsRemoved, sparse (numRuns x numPartitions):\n');
disp(xDimsRemoved_sparse')

%% Frequency domain

fprintf('Convergence, freq (numRuns x numPartitions):\n');
disp(Convergence_freq')
fprintf('DecreasingLowerBound, freq (numRuns x numPartitions):\n');
disp(DecreasingLowerBound_freq')
fprintf('PrivateVarianceFloor, freq (numRuns x numPartitions):\n');
disp(PrivateVarianceFloor_freq')
fprintf('xDimsRemoved, freq (numRuns x numPartitions):\n');
disp(xDimsRemoved_freq')