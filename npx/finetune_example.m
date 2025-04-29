% finetune_example.m
%
% Description: Inspect the effect of using mDLAG-time to fine-tune a model
%              fit by mDLAG-frequency. Only latents found to be
%              significantly shared between areas V1 and V2 are shown.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% Display fitting progress

load('./results/finetune_example.mat');
plotFittingProgress(trackedParams, binWidth, covType, ...
                    'freqLB', freqLB, ...
                    'freqParam', freqParam, ...
                    'units', units);
