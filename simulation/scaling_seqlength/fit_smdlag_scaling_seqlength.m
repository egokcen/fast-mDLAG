% fit_smdlag_scaling_seqlength.m
%
% Description: Fit mDLAG models via inducing variables (smDLAG) to data
%              with increasing sequence (trial) lengths.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_smdlag_scaling_seqlength;

%% Fit smDLAG models

% Set up parallelization
StartParPool(numWorkers);

estParamsList = cell(numRuns,numPartitions);
trackedParamsList = cell(numRuns,numPartitions);
flagsList = cell(numRuns,numPartitions);
parfor runIdx = 1:numRuns
    
    % Load synthetic data and ground truth
    currDataFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
    ws = load(currDataFile);
    
    % Sweep over increasing trial lengths
    estParamsRun = cell(1,numPartitions);
    trackedParamsRun = cell(1,numPartitions);
    flagsRun = cell(1,numPartitions);

    for i = 1:numPartitions
        
        T = Tlist(i);
        fprintf('Run %d of %d, T = %d...\n', runIdx, numRuns, T);

        % Subsample a portion of the trial
        seqSub = [];
        for n = train
            seqSub(n).trialId = ws.seqTrue(n).trialId;
            seqSub(n).T = T;
            seqSub(n).y = ws.seqTrue(n).y(:,1:T);
        end

        % We'll use as few inducing points as are needed to maintain the
        % Nyquist rate for these data (easy to determine, since we know the
        % ground truth timescales).
        S = ceil(T * binWidth / nyqperiod);

        % Initialize model
        initParams = init_smdlag(seqSub, ...
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
                        seqSub, ...
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
        
        T = Tlist(i);
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
        currSaveFile = sprintf('%s/T%04d_sparse.mat', currSaveDir, T);
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end