% pred_summary.m
%
% Description: Compare predictive performance across mDLAG-time,
%              mDLAG-inducing, and mDLAG-frequency.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu

%% mDLAG-time vs mDLAG-frequency and mDLAG-inducing with 20 inducing points

load('./results/pred_mdlag_time.mat');
load('./results/pred_mdlag_freq.mat');
load('./results/pred_smdlag_S20.mat');

figure;
hold on;
line([0 0.1], [0 0.1], 'color', 'k', 'linestyle', '--');
plot(R2mdlag_time(:), R2smdlag(:), ...
    'linestyle', 'none', ...
    'color', SPARSECOLOR, ...
    'marker', '.', ...
    'markersize', 10);
plot(R2mdlag_time(:), R2mdlag_freq(:), ...
    'linestyle', 'none', ...
    'color', FREQCOLOR, ...
    'marker', '.', ...
    'markersize', 10);
axis square;
xlabel('Leave-group-out R^2, mDLAG-time');
ylabel('Leave-group-out R^2, other');

% Significance tests
[p1,~] = signtest(R2mdlag_time(:), R2smdlag(:), ...
    'alpha', 0.05, 'tail', 'right');
[p2,~] = signtest(R2mdlag_freq(:), R2smdlag(:), ...
    'alpha', 0.05, 'tail', 'right');
fprintf('mDLAG-time better than mDLAG-inducing: p = %0.04f\n', p1);
fprintf('mDLAG-freq better than mDLAG-inducing: p = %0.04f\n', p2);

%% mDLAG-time vs mDLAG-frequency and mDLAG-inducing with 32 inducing points

load('./results/pred_mdlag_time.mat');
load('./results/pred_mdlag_freq.mat');
load('./results/pred_smdlag_S32.mat');

figure;
hold on;
line([0 0.1], [0 0.1], 'color', 'k', 'linestyle', '--');
plot(R2mdlag_time(:), R2smdlag(:), ...
    'linestyle', 'none', ...
    'color', SPARSECOLOR, ...
    'marker', '.', ...
    'markersize', 10);
plot(R2mdlag_time(:), R2mdlag_freq(:), ...
    'linestyle', 'none', ...
    'color', FREQCOLOR, ...
    'marker', '.', ...
    'markersize', 10);
axis square;
xlabel('Leave-group-out R^2, mDLAG-time');
ylabel('Leave-group-out R^2, other');

% Significance tests
[p1,~] = signtest(R2mdlag_time(:), R2smdlag(:), ...
    'alpha', 0.05, 'tail', 'right');
[p2,~] = signtest(R2mdlag_freq(:), R2smdlag(:), ...
    'alpha', 0.05, 'tail', 'right');
fprintf('mDLAG-time better than mDLAG-inducing: p = %0.04f\n', p1);
fprintf('mDLAG-freq better than mDLAG-inducing: p = %0.04f\n', p2);
