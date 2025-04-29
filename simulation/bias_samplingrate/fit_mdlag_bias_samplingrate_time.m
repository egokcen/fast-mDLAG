% fit_mdlag_bias_samplingrate_time.m
%
% Description: Fit mDLAG models via time domain to data with increasing 
%              sampling rates.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_bias_samplingrate;

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
    
    % Sweep over increasing sampling rates
    estParamsRun = cell(1,numPartitions);
    trackedParamsRun = cell(1,numPartitions);
    flagsRun = cell(1,numPartitions);

    for i = 1:numPartitions
        
        binWidth = binWidthList(i);
        fprintf('Run %d of %d, sampling period = %d ms...\n', runIdx, numRuns, binWidth);

        % Downsample each trial
        sampleratio = floor(binWidth / min(binWidthList));
        seqSub = [];
        for n = train
            seqSub(n).trialId = ws.seqTrue(n).trialId;
            seqSub(n).y = ws.seqTrue(n).y(:,1:sampleratio:Ttotal);
            seqSub(n).T = size(seqSub(n).y,2);
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
            = em_mdlag(initParams, ...
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
        
        binWidth = binWidthList(i);
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
        currSaveFile = sprintf('%s/period%03d_time.mat', currSaveDir, binWidth);
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end