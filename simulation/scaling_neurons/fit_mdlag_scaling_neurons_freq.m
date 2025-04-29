% fit_mdlag_scaling_neurons_freq.m
%
% Description: Fit mDLAG models via frequency domain approximation to 
%              data with increasing neuron number.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_scaling_neurons;

%% Fit mDLAG models

% Set up parallelization
StartParPool(numWorkers);

estParamsList = cell(numRuns,numPartitions);
trackedParamsList = cell(numRuns,numPartitions);
flagsList = cell(numRuns,numPartitions);
parfor runIdx = 1:numRuns
    
    % Load synthetic data and ground truth
    currDataFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
    ws = load(currDataFile);
    
    % Sweep over increasing neuron number
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

        % Initialize model
        initParams = init_mdlag(seqSub, ...
                                yDims, ...
                                xDim_fit, ...
                                binWidth,...
                                'covType', covType, ...
                                'prior', prior, ...
                                'randomSeed', randomSeed, ...
                                'saveCcov', saveCcov);
        
        % As a preprocessing step, compute the FFT of the observed data
        seqSub = fftseq(seqSub, 'y', 'yfft');
        
        % Fit model
        [estParams,~,trackedParams,flags] ...
            = em_mdlag_freq(initParams, ...
                            seqSub, ...
                            xDim_fit, ...
                            'prior', prior, ...
                            'tol', tol, ...
                            'maxIters', maxIters, ...
                            'freqLB', freqLB, ...
                            'freqParam', freqParam, ...
                            'learnDelays', learnDelays, ...
                            'verbose', verbose, ...
                            'minVarFrac', minVarFrac, ...
                            'maxDelayFrac', maxDelayFrac, ...
                            'maxTauFrac', maxTauFrac, ...
                            'pruneX', pruneX, ...
                            'saveXcov', saveXcov, ...
                            'saveCcov', saveCcov);
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
        currSaveFile = sprintf('%s/yDim%03d_freq.mat', currSaveDir, sum(yDims));
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end