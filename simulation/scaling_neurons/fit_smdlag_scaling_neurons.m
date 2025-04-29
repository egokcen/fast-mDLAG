% fit_smdlag_scaling_neurons.m
%
% Description: Fit mDLAG models via inducing variables (smDLAG) to data
%              with increasing neuron number.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_smdlag_scaling_neurons;

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
        
        yDims = yDimList(i,:);
        fprintf('Run %d of %d, yDim = %d...\n', runIdx, numRuns, sum(yDims));

        % Take observations from only a subset of neurons
        seqSub = [];
        seqSub_grouped = partitionSeq(ws.seqTrue, yDimList(end,:), 'datafield', 'y');
        for n = train
            seqSub(n).y = [];
            seqSub(n).trialId = ws.seqTrue(n).trialId;
            seqSub(n).T = ws.seqTrue(n).T;
            for groupIdx = 1:numGroups
                seqSub(n).y = [seqSub(n).y; seqSub_grouped{groupIdx}(n).y(1:yDims(groupIdx),:)];
            end
        end

        % We'll use as few inducing points as are needed to maintain the
        % Nyquist rate for these data (easy to determine, since we know the
        % ground truth timescales).
        S = ceil(Ttotal * binWidth / nyqperiod);

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
        
        yDims = yDimList(i,:);
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
        currSaveFile = sprintf('%s/yDim%03d_sparse.mat', currSaveDir, sum(yDims));
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end