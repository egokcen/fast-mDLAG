% generate_data_bias_samplingrate.m
%
% Description: Generate synthetic datasets to characterize parameter
%              estimation bias as a function of sampling rate. Ground 
%              truth data and parameters are saved to files.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_bias_samplingrate;

%% Randomly generate data from mDLAG model

rng(0); % Seed the random number generator for reproducibility
for runIdx = 1:numRuns

    fprintf('Run %d of %d\n', runIdx, numRuns);

    % Generate data at the highest sampling rate, and then during analyses
    % we'll downsample the data.

    % Even data generation in the time domain is prohibitively expensive.
    % Conservatively generate 3x as much data via the frequency domain, and
    % then only take the middle Ttotal time points of the trial to 
    % effectively remove periodic edge effects.
    [seqTrue, trueParams] ...
        = simdata_mdlag_freq(Ntotal, 3*Ttotal, min(binWidthList), yDims, xDim, ...
                             hyperparams, snr, covType, gp_params);
    % Remove periodic edge effects
    seqTrue = rmfield(seqTrue, 'xfft');
    for n = 1:length(seqTrue)
        seqTrue(n).T = Ttotal;
        seqTrue(n).xsm = seqTrue(n).xsm(:,Ttotal+1:2*Ttotal);
        seqTrue(n).y = seqTrue(n).y(:,Ttotal+1:2*Ttotal);
    end
    
    % Save generated data, along with ground truth parameters
    currSaveDir = sprintf('%s/run%03d', dataDir, runIdx);
    if isfolder(currSaveDir)
        fprintf('Using existing directory %s...\n', currSaveDir);
    else
        fprintf('Making directory %s...\n', currSaveDir);
        mkdir(currSaveDir);
    end
    currSaveFile = sprintf('%s/%s', currSaveDir, dataFile);
    save(currSaveFile, 'seqTrue', 'trueParams');
end

%% Inspect generated data parameters

clear jointParams;

for runIdx = 1:numRuns
    
    % Load generated data, along with ground truth parameters
    currSaveFile = sprintf('%s/run%03d/%s', dataDir, runIdx, dataFile);
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

plotGPparams_mdlag(jointParams,min(binWidthList),...
                   'sigDims',repmat(sigDimsTrue,1,numRuns),...
                   'units',units);
