%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpaRWE: Reliability Quantification and Weighted Ensemble
% Single-slice demonstration
%
% Input:
%   data/507.mat, containing a variable named "baseCls"
%
% Output:
%   results/SpaRWE_507_all_betas.csv
%
% The algorithmic calculations are retained from the supplied source code.
% Function names, file organization, and the demonstration scope are
% simplified for the SpaRWE code repository.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ===============================
% User settings
%% ===============================

beta_list = [0, 1, 2, 4, 8];

% Demonstration slice and its number of clusters
slice_id = 507;
clsNums = 7;

% Use paths relative to this script
project_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(project_dir, 'data');
results_dir = fullfile(project_dir, 'results');

mat_file = fullfile(data_dir, sprintf('%d.mat', slice_id));

if exist(mat_file, 'file') ~= 2
    error('Input file not found: %s\nPlace %d.mat in the data folder.', ...
        mat_file, slice_id);
end

loaded_data = load(mat_file, 'baseCls');

if ~isfield(loaded_data, 'baseCls')
    error('The input MAT-file must contain a variable named "baseCls".');
end

baseCls = double(loaded_data.baseCls);

fprintf('========================================\n');
fprintf('SpaRWE: processing slice %d\n', slice_id);
fprintf('Samples: N = %d\n', size(baseCls, 1));
fprintf('Base clusterings: M = %d\n', size(baseCls, 2));
fprintf('Number of clusters: K = %d\n', clsNums);
fprintf('========================================\n');

%% ===============================
% Part 1. Reliability quantification
%% ===============================

fprintf('\nCalculating reliability weights...\n');

ncai = calculate_reliability(baseCls);

%% ===============================
% Part 2. Weighted ensemble
%% ===============================

fprintf('\nRunning weighted ensemble...\n');

all_results = [];

for b_idx = 1:length(beta_list)

    beta = beta_list(b_idx);
    fprintf('  beta = %g\n', beta);

    I_ncai = ncai.^beta;

    % Use the Average Linkage result (third output)
    resultCls = weighted_ensemble(baseCls, I_ncai, clsNums);

    all_results = [all_results, resultCls];

end

%% ===============================
% Save final clustering results
%% ===============================

if exist(results_dir, 'dir') ~= 7
    mkdir(results_dir);
end

beta_names = cell(1, length(beta_list));

for b_idx = 1:length(beta_list)
    beta_names{b_idx} = sprintf('beta_%g', beta_list(b_idx));
end

results_table = array2table(all_results, ...
    'VariableNames', beta_names);

output_file = fullfile(results_dir, ...
    sprintf('SpaRWE_%d_all_betas.csv', slice_id));

writetable(results_table, output_file);

fprintf('\nSaved: %s\n', output_file);
fprintf('SpaRWE completed.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1. Reliability quantification functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NCAI = calculate_reliability(baseCls)

cnt = size(baseCls, 2);
ClsAfnyWrtNMI = zeros(cnt, cnt); % It is symmetric.

for i = 1:cnt
    for j = i+1:cnt
        ClsAfnyWrtNMI(i,j) = nmi_max(baseCls(:,i), baseCls(:,j));
    end
end

ClsAfnyWrtNMI = ClsAfnyWrtNMI + ClsAfnyWrtNMI';

CAI = (sum(ClsAfnyWrtNMI)./(cnt-1))';
NCAI = CAI./max(CAI(:));

end


function value = nmi_max(x, y)

% Compute normalized mutual information: I(x,y) / max(H(x), H(y)).

assert(numel(x) == numel(y));

n = numel(x);
x = reshape(x, 1, n);
y = reshape(y, 1, n);

l = min(min(x), min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x), max(y));

idx = 1:n;
Mx = sparse(idx, x, 1, n, k, n);
My = sparse(idx, y, 1, n, k, n);
Pxy = nonzeros(Mx'*My/n); % Joint distribution of x and y
Hxy = -dot(Pxy, log(Pxy+eps));

Px = mean(Mx, 1);
Py = mean(My, 1);

% Entropy of Py and Px
Hx = -dot(Px, log(Px+eps));
Hy = -dot(Py, log(Py+eps));

% Mutual information
MI = Hx + Hy - Hxy;

% Maximum normalized mutual information
value = MI/max(Hx, Hy);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2. Weighted ensemble functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resultCls = weighted_ensemble(baseCls, I_ncai, clsNums)

% Get the weighted co-association matrix.
S = get_weighted_matrix(baseCls, I_ncai);

d = similarity_to_distance(S);
clear S

% Average Linkage
Zal = linkage(d, 'average');
clear d

resultCls = [];

for K = clsNums
    resultCls = [resultCls cluster(Zal, 'maxclust', K)];
end

end


function S = get_weighted_matrix(baseCls, I_ncai)

I_ncai = I_ncai(:);
w_sum = sum(I_ncai);

[n, nBC] = size(baseCls);
cntCol = max(baseCls);

S = zeros(n, n);

for k = 1:nBC
    tmp = double(bsxfun(@eq, 1:cntCol(k), baseCls(:,k)));
    idx = tmp * tmp';
    S = S + idx.*I_ncai(k);
end

S = S/w_sum;

end


function d = similarity_to_distance(S)

% Convert the similarity matrix to a distance vector for linkage.

s = [];

for a = 1:length(S)-1
    s = [s S(a, [a+1:end])];
end

d = 1 - s;

end
