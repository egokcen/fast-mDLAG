% fit_mdlag_scaling_numgroups_freq_taper.m
%
% Description: Fit mDLAG models via frequency domain approximation to 
%              data with increasing group number. Apply a taper to the 
%              training data to see if it mitigates bias in GP parameter 
%              estimates.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_scaling_numgroups;

%% Fit mDLAG models

% Set up parallelization
StartParPool(numWorkers);

estParamsList = cell(numRuns,numPartitions);
trackedParamsList = cell(numRuns,numPartitions);
flagsList = cell(numRuns,numPartitions);
parfor runIdx = 1:numRuns
    
    % Sweep over increasing trial lengths
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

        % Apply the taper
        seqTrain = apply_taper(seqTrain, @(T) hamming(T, 'periodic'), 'y');

        % Initialize model
        initParams = init_mdlag(seqTrain, ...
                                yDims, ...
                                xDim_fit, ...
                                binWidth,...
                                'covType', covType, ...
                                'prior', prior, ...
                                'randomSeed', randomSeed, ...
                                'saveCcov', saveCcov);
        
        % As a preprocessing step, compute the FFT of the observed data
        seqTrain = fftseq(seqTrain, 'y', 'yfft');
        
        % Fit model
        [estParams,~,trackedParams,flags] ...
            = em_mdlag_freq(initParams, ...
                            seqTrain, ...
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
        currSaveFile = sprintf('%s/numgroups_%02d_freq_hamming.mat', currSaveDir, numGroups);
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end