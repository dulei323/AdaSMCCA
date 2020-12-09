
close all; clear; clc;
% -----------------------------------
% Author: Lei Du, dulei@nwpu.edu.cn
% Date: 09-Dec-2020
% -----------------------------------

%% Load data
load('synData.mat');
for i = 1 : numel(Data_X)
    X{i} = Data_X{i};
end

%% Tuned parameters
% Robustness-aware AdaSMCCA
opts.rAdaSMCCA.lambda_1 = 1;
opts.rAdaSMCCA.beta = 0.1;    % FGL & L1-norm for SNP
opts.rAdaSMCCA.lambda_2 = 1;  % L1-norm for Protein
opts.rAdaSMCCA.lambda_3 = 1;  % L1-norm for Imaging
% Uncertainty-aware AdaSMCCA
opts.unAdaSMCCA.lambda_1 = 100;
opts.unAdaSMCCA.beta = 0.1;      % FGL & L1-norm for SNP
opts.unAdaSMCCA.lambda_2 = 100;  % L1-norm for Protein
opts.unAdaSMCCA.lambda_3 = 100;  % L1-norm for Imaging

%% Kfold Cross validation
n = size(X{1}, 1);
k_fold = 5;
indices = crossvalind('Kfold', n, k_fold);

fprintf('===================================\n');
for k = 1 : k_fold
    fprintf('Current fold: %d\n', k);
    
    % Split training data and test data
    test = (indices == k);
    train = ~test;
    for i = 1 : numel(X)
        trainData.X{i} = normalize(X{i}(train, :), 'norm');
        testData.X{i} = normalize(X{i}(test, :), 'norm');
    end
    
    % Training step
    % Robustness-aware AdaSMCCA
    tic;
    [W.rAdaSMCCA{k}, u.rAdaSMCCA(:, k), v.rAdaSMCCA(:, k), w.rAdaSMCCA(:, k)] = rAdaSMCCA(trainData, opts.rAdaSMCCA);
    fprintf('Robustness-aware AdaSMCCA: %.3fs\n', toc);
    % Uncertainty-aware AdaSMCCA
    tic;
    [W.unAdaSMCCA{k}, u.unAdaSMCCA(:, k), v.unAdaSMCCA(:, k), w.unAdaSMCCA(:, k)] = unAdaSMCCA(trainData, opts.unAdaSMCCA);
    fprintf('Uncertainty-aware AdaSMCCA: %.3fs\n', toc);
    
    % Canonical Correlation Coefficients (CCCs)
    % Robustness-aware AdaSMCCA
    CCC_train.rAdaSMCCA(k, :) = calcCCC(trainData, W.rAdaSMCCA{k});
    CCC_test.rAdaSMCCA(k, :) = calcCCC(testData, W.rAdaSMCCA{k});
    % Uncertainty-aware AdaSMCCA
    CCC_train.unAdaSMCCA(k, :) = calcCCC(trainData, W.unAdaSMCCA{k});
    CCC_test.unAdaSMCCA(k, :) = calcCCC(testData, W.unAdaSMCCA{k});
    
    if k ~= k_fold
        fprintf('\n');
    end
end
fprintf('===================================\n');

%% Weights
% Robustness-aware AdaSMCCA
u.rAdaSMCCA_mean = mean(u.rAdaSMCCA, 2);
v.rAdaSMCCA_mean = mean(v.rAdaSMCCA, 2);
w.rAdaSMCCA_mean = mean(w.rAdaSMCCA, 2);
% Uncertainty-aware AdaSMCCA
u.unAdaSMCCA_mean = mean(u.unAdaSMCCA, 2);
v.unAdaSMCCA_mean = mean(v.unAdaSMCCA, 2);
w.unAdaSMCCA_mean = mean(w.unAdaSMCCA, 2);

%% CCCs
% Robustness-aware AdaSMCCA
CCC_train.rAdaSMCCA_mean = mean(CCC_train.rAdaSMCCA, 1);
CCC_test.rAdaSMCCA_mean = mean(CCC_test.rAdaSMCCA, 1);
% Uncertainty-aware AdaSMCCA
CCC_train.unAdaSMCCA_mean = mean(CCC_train.unAdaSMCCA, 1);
CCC_test.unAdaSMCCA_mean = mean(CCC_test.unAdaSMCCA, 1);

%% Figures
figure; colormap('jet');

% Ground Truth
subplot(3, 3, 1); imagesc(GroundTruth_W{1}', [-3 3]);
set(gca, 'XTick', [], 'YTick', []);
ylabel('Ground Truth');
title('w1');
subplot(3, 3, 2); imagesc(GroundTruth_W{2}', [-3 3]);
set(gca, 'XTick', [], 'YTick', []);
title('w2');
subplot(3, 3, 3); imagesc(GroundTruth_W{3}', [-3 3]);
set(gca, 'XTick', [], 'YTick', []);
title('w3');
colorbar('Ticks', [-3 0 3]);

% Robustness-aware AdaSMCCA
subplot(3, 3, 4); imagesc(u.rAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
ylabel('rAdaSMCCA');
subplot(3, 3, 5); imagesc(v.rAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
subplot(3, 3, 6); imagesc(w.rAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
colorbar('Ticks', [-0.1 0 0.1]);

% Uncertainty-aware AdaSMCCA
subplot(3, 3, 7); imagesc(u.unAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
ylabel('unAdaSMCCA');
subplot(3, 3, 8); imagesc(v.unAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
subplot(3, 3, 9); imagesc(w.unAdaSMCCA_mean', [-0.1 0.1]);
set(gca, 'XTick', [], 'YTick', []);
colorbar('Ticks', [-0.1 0 0.1]);
