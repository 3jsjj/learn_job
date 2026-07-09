clear;
clc;
close all;

%% 1. 加载数据库
if ~isfile('bayesian_inverse_database.mat')
    error('请先运行 build_bayesian_database.m。');
end

load('bayesian_inverse_database.mat');

%% 2. 输入新的9个传感器值
% 顺序：
% S1=(0,0), S2=(10,0), S3=(20,0)
% S4=(0,10),S5=(10,10),S6=(20,10)
% S7=(0,20),S8=(10,20),S9=(20,20)

observed_sensor = [
    0;
    0;
    0;
    0;
    1000;
    0;
    0;
    0;
    0
];

observed_sensor = observed_sensor(:)';

if numel(observed_sensor) ~= 9
    error('observed_sensor 必须包含9个值。');
end

%% 3. 贝叶斯反演
result = bayesian_inverse_predict( ...
    observed_sensor, ...
    X_sensor, ...
    Y_pressure, ...
    sigma_sensor, ...
    prior_probability);

%% 4. MAP结果
idx_map = result.map_index;

fprintf('\n================ 概率反演结果 ================\n');
fprintf('MAP工况：%s\n', case_name(idx_map));
fprintf('最可能压头边长：%.4f mm\n', side_length_mm(idx_map));
fprintf('最可能压头面积：%.4f mm^2\n', area_mm2(idx_map));
fprintf('最可能压入深度：%.4f mm\n', depth_mm(idx_map));
fprintf('MAP后验概率：%.8f\n', result.posterior_probability(idx_map));
fprintf('后验熵：%.6f\n', result.posterior_entropy);
fprintf('有效候选数：%.4f\n', result.effective_sample_size);
fprintf('==============================================\n');

%% 5. 后验均值参数
posterior_side = sum(result.posterior_probability .* side_length_mm);
posterior_area = sum(result.posterior_probability .* area_mm2);
posterior_depth = sum(result.posterior_probability .* depth_mm);

fprintf('后验均值边长：%.4f mm\n', posterior_side);
fprintf('后验均值面积：%.4f mm^2\n', posterior_area);
fprintf('后验均值深度：%.4f mm\n\n', posterior_depth);

%% 6. 前10个候选工况
top_k = min(10, numel(result.posterior_probability));
[sorted_prob, sorted_idx] = sort(result.posterior_probability, 'descend');

TopTable = table( ...
    case_name(sorted_idx(1:top_k)), ...
    side_length_mm(sorted_idx(1:top_k)), ...
    area_mm2(sorted_idx(1:top_k)), ...
    depth_mm(sorted_idx(1:top_k)), ...
    sorted_prob(1:top_k), ...
    'VariableNames', { ...
    'CaseName', ...
    'SideLength_mm', ...
    'Area_mm2', ...
    'Depth_mm', ...
    'PosteriorProbability'});

disp(TopTable);
writetable(TopTable, 'bayesian_top_candidates.csv');

%% 7. 压力场
F_map_map = reshape(result.map_pressure, 21, 21);
F_map_mean = reshape(result.posterior_mean_pressure, 21, 21);
F_map_std = reshape(result.posterior_std_pressure, 21, 21);

ResultTable = table( ...
    grid_points(:,1), ...
    grid_points(:,2), ...
    result.map_pressure(:), ...
    result.posterior_mean_pressure(:), ...
    result.posterior_std_pressure(:), ...
    'VariableNames', { ...
    'x_mm', ...
    'y_mm', ...
    'MAP_CPRESS', ...
    'PosteriorMean_CPRESS', ...
    'PosteriorStd_CPRESS'});

writetable(ResultTable, 'bayesian_pressure_reconstruction.csv');

%% 8. 绘图
figure('Name','Bayesian MAP Pressure Field');
surf(X_grid, Y_grid, F_map_map, 'EdgeColor','none');
colormap jet;
colorbar;
xlabel('x / mm');
ylabel('y / mm');
zlabel('CPRESS');
title('Bayesian MAP reconstructed pressure field');
axis tight;
grid on;
view(45,30);

figure('Name','Bayesian Posterior Mean Pressure Field');
surf(X_grid, Y_grid, F_map_mean, 'EdgeColor','none');
colormap jet;
colorbar;
xlabel('x / mm');
ylabel('y / mm');
zlabel('CPRESS');
title('Bayesian posterior mean pressure field');
axis tight;
grid on;
view(45,30);

figure('Name','Bayesian Uncertainty');
contourf(X_grid, Y_grid, F_map_std, 30, 'LineColor','none');
colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
axis tight;
title('Posterior standard deviation');

figure('Name','Posterior Probability');
bar(sorted_prob(1:top_k));
xlabel('Candidate rank');
ylabel('Posterior probability');
title('Top posterior probabilities');
grid on;
