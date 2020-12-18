clc; close all; clear;

load('adj.mat');

%% Collect edge weights

% In the WSBM code missing values are represented by NaNs.
% However, functions of the BCT typically need 0s for missing values.
% Hence, conversions of NaNs to 0s are made here.
adj(adj == 0) = NaN;
% Are drugs trivially co-prescribed with themselves?
adj(1:size(adj,1)+1:end) = NaN; % remove diagonal
edge_list = Adj2Edg(adj);
edge_weights = edge_list(:,3);

figure;
hist = histfit(log(edge_weights));
hist(1).FaceColor = [.8 .8 1];
hist(2).Color = [.2 .2 .2];
xlabel('Log Edge Weight', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);

%% Run WSBM

% input parameters for WSBM
W_Distr = 'Normal';
E_Distr = 'Bernoulli';
num_iters = 25;

log_edge_list = edge_list;
log_edge_list(:,3) = log(edge_list(:,3));

for k = 2:1:20 % sweep over number of communities
    ModelInputs = cell(numel(num_iters),1);
    for iter = 1:num_iters
        ModelInputs{iter} = {k, 'W_Distr', W_Distr, 'E_Distr', E_Distr};
    end
    [Best_Model,Scores,Models] = wsbmLooper(log_edge_list, ModelInputs);
    filename = sprintf('k_%d_communities', k);
    save(fullfile('WSBM_results', filename), 'Best_Model', 'Scores', 'Models');
end

%% Determine k

path_1 = '/Users/Shubhankar/Desktop/Developer/CIS-545/WSBM_results';
k = 2:1:20;
log_evidence = zeros(1,length(k));
errors = zeros(1,length(k));

for idx = 1:length(k)
    curr_k = k(idx);
    path_2 = sprintf('k_%d_communities.mat', curr_k);
    load(fullfile(path_1, path_2)); % loads in Scores and some other variables
    log_evidence(idx) = mean(Scores);
    errors(idx) = std(Scores);
end

num_iters = 25;

figure;
e = errorbar(k, log_evidence, errors/sqrt(num_iters), ...
    '.k', 'MarkerSize', 25, 'LineWidth', 0.1,...
    'MarkerEdgeColor','black','MarkerFaceColor','black');
xlabel('Number of Blocks', 'FontSize', 15);
ylabel('Log Likelihood', 'FontSize', 15);
title('Number of Communities vs. Log Likelihood', 'FontSize', 15);
prettify;

%% Get consensus partition

% Optimal number of communities is k = 4
optimal_k = 4;
path_1 = '/Users/Shubhankar/Desktop/Developer/CIS-545/WSBM_results';
path_2 = sprintf('k_%d_communities.mat', optimal_k);
load(fullfile(path_1, path_2)); % loads in Models and some other variables

VI_mat = zeros(length(Models), length(Models));
for i = 1:length(Models)
    model_1 = Models{i, 1}.Para.mu;
    for j = 1:length(Models)
        model_2 = Models{j, 1}.Para.mu;
        VI_mat(i, j) = -varInfo(model_1, model_2);
    end
end
[~, best_model_idx] = max(sum(VI_mat, 2)); % index of most central model
[~, partition] = max(Models{best_model_idx, 1}.Para.mu);

community_1 = find(partition == 1);
community_2 = find(partition == 2);
community_3 = find(partition == 3);
community_4 = find(partition == 4);

[X, Y, INDSORT] = grid_communities(partition); 
figure;
imagesc(adj(INDSORT, INDSORT)); % plot ordered adjacency matrix
hold on; % hold on to overlay community visualization
plot(X, Y, 'r', 'linewidth', 2);  
hold off;



