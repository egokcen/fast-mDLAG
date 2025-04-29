% generate_data_longtrials.m
%
% Description: Generate synthetic datasets with very long trials, to
%              further demonstrate the scalability of mDLAG-frequency.
%              Ground truth data and parameters are saved to files.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_scaling_seqlength;

%% Randomly generate data from mDLAG model

rng(0); % Seed the random number generator for reproducibility
for runIdx = 1:numRuns

    fprintf('Run %d of %d\n', runIdx, numRuns);

    % Even data generation in the time domain is prohibitively expensive.
    % Conservatively generate 20% longer trials via the frequency domain,
    % and then only take the middle Tlong time points of the trial to 
    % effectively remove periodic edge effects.
    [seqTrue, trueParams] ...
        = simdata_mdlag_freq(Ntotal, ceil(1.2*Tlong), binWidth, yDims, xDim, ...
                             hyperparams, snr, covType, gp_params);
    % Remove periodic edge effects
    seqTrue = rmfield(seqTrue, 'xfft');
    for n = 1:length(seqTrue)
        seqTrue(n).T = Tlong;
        seqTrue(n).xsm = seqTrue(n).xsm(:,ceil(0.1*Tlong)+(1:Tlong));
        seqTrue(n).y = seqTrue(n).y(:,ceil(0.1*Tlong)+(1:Tlong));
    end
    
    % Save generated data, along with ground truth parameters
    currSaveDir = sprintf('%s/run%03d', dataDir, runIdx);
    if isfolder(currSaveDir)
        fprintf('Using existing directory %s...\n', currSaveDir);
    else
        fprintf('Making directory %s...\n', currSaveDir);
        mkdir(currSaveDir);
    end
    currSaveFile = sprintf('%s/%s', currSaveDir, longDataFile);
    save(currSaveFile, 'seqTrue', 'trueParams', 'snr', 'binWidth');
end

%% Inspect generated data parameters

clear jointParams;

for runIdx = 1:numRuns
    
    % Load generated data, along with ground truth parameters
    currSaveFile = sprintf('%s/run%03d/%s', dataDir, runIdx, longDataFile);
    load(currSaveFile);
    
    % Collect all GP parameters into one structure so we can plot them
    if runIdx == 1
        % Initialize the joint parameter structures
        jointParams = trueParams;
    else
        % Otherwise, concatenate parameters across datasets
        jointParams.xDim = jointParams.xDim + trueParams.xDim;
        jointParams.D = [jointParams.D trueParams.D];
        jointParams.gamma = [jointParams.gamma trueParams.gamma];
        jointParams.eps = [jointParams.eps trueParams.eps];
    end
 
end

%% Visualize

plotGPparams_mdlag(jointParams,binWidth,...
                   'sigDims',repmat(sigDimsTrue,1,numRuns),...
                   'units',units);
