
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpatialESD: Ensemble-based spatial domain identification
%
% Input: Five base clustering results
% Output: Final spatial domain assignments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ========================= Parameters ==================================

slice_id = '507';

spot_number = [];       % Automatically determined
method_number = 5;      % Number of base clustering methods
K = 7;                  % Number of final spatial domains
t = 20;                 % Random-walk iteration parameter

input_dir = 'D:/your_project/ComparativeMethods/input';
output_dir = 'D:/your_project/ComparativeMethods/results';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ========================= Load base clustering results =================

data1 = read_base_labels(fullfile(input_dir, ['data1_' slice_id '.csv']));
data2 = read_base_labels(fullfile(input_dir, ['data2_' slice_id '.csv']));
data3 = read_base_labels(fullfile(input_dir, ['data3_' slice_id '.csv']));
data4 = read_base_labels(fullfile(input_dir, ['data4_' slice_id '.csv']));
data5 = read_base_labels(fullfile(input_dir, ['data5_' slice_id '.csv']));

%% ========================= Construct baseCls ============================

spot_number = length(data1);

if any([length(data2), length(data3), ...
        length(data4), length(data5)] ~= spot_number)
    error('The five base clustering results have different spot numbers.');
end

baseCls = zeros(spot_number, method_number);

baseCls(:,1) = data1;
baseCls(:,2) = data2;
baseCls(:,3) = data3;
baseCls(:,4) = data4;
baseCls(:,5) = data5;

% Convert cluster labels into consecutive positive integers
for j = 1:method_number
    [~,~,baseCls(:,j)] = unique(baseCls(:,j), 'stable');
end

%% ========================= SpatialESD ensemble ==========================

ensemble = SpatialESD(baseCls, K, t);

%% ========================= Save results =================================

result = table( ...
    (1:spot_number)', ...
    ensemble(:), ...
    'VariableNames', {'spot_index', 'spatial_domain'});

output_file = fullfile(output_dir, ...
    ['SpatialESD_' slice_id '.csv']);

writetable(result, output_file);

fprintf('SpatialESD completed!\n');
fprintf('Number of spots: %d\n', spot_number);
fprintf('Number of base methods: %d\n', method_number);
fprintf('Number of spatial domains: %d\n', ...
    length(unique(ensemble)));
fprintf('Results saved to: %s\n', output_file);


%% ========================= Data loading function ========================

function labels = read_base_labels(file)

    if ~isfile(file)
        error('File not found: %s', file);
    end

    T = readtable(file, 'VariableNamingRule', 'preserve');

    labels = T{:, end};
    labels = double(labels(:));

    if isempty(labels) || any(~isfinite(labels))
        error('Invalid cluster labels in: %s', file);
    end

end
