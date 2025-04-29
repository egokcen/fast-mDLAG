% fit_smdlag_scaling_numgroups.m
%
% Description: Fit mDLAG models via inducing variables (smDLAG) to data
%              with increasing group number.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_smdlag_scaling_numgroups;

%% Fit smDLAG models

% Set up parallelization
StartParPool(numWorkers);

estParamsList = cell(numRuns,numPartitions);
trackedParamsList = cell(numRuns,numPartitions);
flagsList = cell(numRuns,numPartitions);
parfor runIdx = 1:numRuns
    
    % Sweep over increasing group numbers
    estParamsRun = cell(1,numPartitions);
    trackedParamsRun = cell(1,numPartitions);
    flagsRun = cell(1,numPartitions);

    for i = 1:numPartitions
        
        numGroups = numGroupsList(i);
        fprintf('Run %d of %d, M = %d...\n', runIdx, numRuns, numGroups);

        % Load synthetic data and ground truth
        currDataFile = sprintf('%s/run%03d/%s%02d', dataDir, runIdx, dataFile, numGroups);
        ws = load(currDataFile);
        seqTrain = ws.seqTrue(train);
        yDims = ws.trueParams.yDims;

        % Initialize model
        initParams = init_smdlag(seqTrain, ...
                                 yDims, ...
                                 xDim_fit, ...
                                 binWidth,...
                                 'S', S, ...
                                 'prior', prior, ...
                                 'randomSeed', randomSeed, ...
                                 'saveCcov', saveCcov);
        
        % Fit model
        [estParams,~,trackedParams,flags] ...
            = em_smdlag(initParams, ...
                        seqTrain, ...
                        xDim_fit, ...
                        'prior', prior, ...
                        'tol', tol, ...
                        'maxIters', maxIters, ...
                        'freqLB', freqLB, ...
                        'freqParam', freqParam, ...
                        'learnTau', learnTau, ...
                        'learnDelays', learnDelays, ...
                        'learnInducingLocs', learnInducingLocs, ...
                        'verbose', verbose, ...
                        'minVarFrac', minVarFrac, ...
                        'maxDelayFrac', maxDelayFrac, ...
                        'maxTauFrac', maxTauFrac, ...
                        'pruneX', pruneX, ...
                        'saveCcov', saveCcov, ...
                        'saveWcov', saveWcov);
         estParamsRun{i} = estParams;
         trackedParamsRun{i} = trackedParams;
         flagsRun{i} = flags;
                            
    end
    estParamsList(runIdx,:) = estParamsRun;
    trackedParamsList(runIdx,:) = trackedParamsRun;
    flagsList(runIdx,:) = flagsRun;
end

% Save results
for runIdx = 1:numRuns
    
    for i = 1:numPartitions
        
        numGroups = numGroupsList(i);
        estParams = estParamsList{runIdx,i};
        trackedParams = trackedParamsList{runIdx,i};
        flags = flagsList{runIdx,i};
        
        currSaveDir = sprintf('%s/run%03d', resultDir, runIdx);
        if isfolder(currSaveDir)
            fprintf('Using existing directory %s...\n', currSaveDir);
        else
            fprintf('Making directory %s...\n', currSaveDir);
            mkdir(currSaveDir);
        end
        currSaveFile = sprintf('%s/numgroups_%02d_sparse.mat', currSaveDir, numGroups);
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end