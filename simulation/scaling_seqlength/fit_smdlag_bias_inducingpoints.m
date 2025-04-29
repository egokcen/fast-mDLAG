% fit_smdlag_bias_inducingpoints.m
%
% Description: Fit smDLAG models to data using an increasing number of 
%              inducing points.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_smdlag_scaling_seqlength;

%% Fit smDLAG models

% Set up parallelization
StartParPool(numWorkers);

estParamsList = cell(numRuns,numS);
trackedParamsList = cell(numRuns,numS);
flagsList = cell(numRuns,numS);
parfor runIdx = 1:numRuns
    
    % Load synthetic data and ground truth
    currDataFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
    ws = load(currDataFile);
    
    % Sweep over increasing number of inducing points
    estParamsRun = cell(1,numS);
    trackedParamsRun = cell(1,numS);
    flagsRun = cell(1,numS);

    for i = 1:numS
        
        S = Slist(i);
        fprintf('Run %d of %d, S = %d...\n', runIdx, numRuns, S);

        % Subsample a portion of the trial
        seqSub = [];
        for n = train
            seqSub(n).trialId = ws.seqTrue(n).trialId;
            seqSub(n).T = Tfix;
            seqSub(n).y = ws.seqTrue(n).y(:,1:Tfix);
        end

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
    
    for i = 1:numS
        
        S = Slist(i);
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
        currSaveFile = sprintf('%s/T%04d_sparse_S%03d.mat', currSaveDir, Tfix, S);
        save(currSaveFile, 'estParams', 'trackedParams', 'flags');
         
    end
    
end