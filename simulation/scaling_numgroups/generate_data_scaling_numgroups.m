% generate_data_scaling_numgroups.m
%
% Description: Generate synthetic datasets to characterize runtime and 
%              performance as a function of group number. Ground 
%              truth data and parameters are saved to files.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Set up script parameters

set_consts_mdlag_scaling_numgroups;

%% Randomly generate data from mDLAG model

rng(1); % Seed the random number generator for reproducibility
for runIdx = 1:numRuns

    fprintf('Run %d of %d\n', runIdx, numRuns);

    for i = 1:numPartitions

        fprintf('    Partition %d of %d\n', i, numPartitions);
        numGroups = numGroupsList(i);

        % Determine how many units are in each group
        yDims = repmat(floor(yDim ./ numGroups), 1, numGroups);
        block_idxs = get_block_idxs(yDims);

        % ARD parameters
        hyperparams.a_alpha = MAG.*ones(numGroups,xDim);
        hyperparams.b_alpha = MAG.*ones(numGroups,xDim);

        % Signal-to-noise ratios of each group
        curr_snr = snr .* ones(1,numGroups);

        % GP parameters
        gp_params = generate_GPparams_mdlag(xDim, numGroups, covType, lims);

        % Even data generation in the time domain is prohibitively expensive.
        % Conservatively generate 3x as much data via the frequency domain, and
        % then only take the middle Ttotal time points of the trial to 
        % effectively remove periodic edge effects.
        [seqTrue, trueParams] ...
        = simdata_mdlag_freq(Ntotal, 3*Ttotal, binWidth, yDims, xDim, ...
                             hyperparams, curr_snr, covType, gp_params);

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
            fprintf('        Using existing directory %s...\n', currSaveDir);
        else
            fprintf('        Making directory %s...\n', currSaveDir);
            mkdir(currSaveDir);
        end
        currSaveFile = sprintf('%s/%s%02d.mat', currSaveDir, dataFile, numGroups);
        save(currSaveFile, 'seqTrue', 'trueParams', 'binWidth');

    end

end
