clear;
clc;
close all;

rng(0);

%% ============================================================
%  1. 定义传感器位置和当前测量值
% ============================================================

sensor_points = [
    0, 0;
    10, 0;
    20, 0;
    0, 10;
    10, 10;
    20, 10;
    0, 20;
    10, 20;
    20, 20
];

current_sensor_values = [
    1085.53;
    2357.69;
    1085.53;
    2357.69;
    6258.71;
    2357.69;
    1085.53;
    2357.69;
    1085.53
];

%% ============================================================
%  2. 构建 21 × 21 插值网格
% ============================================================

x_int = 0:1:20;
y_int = 0:1:20;

[X_grid, Y_grid] = meshgrid(x_int, y_int);

nx = length(x_int);
ny = length(y_int);
n_grid = nx * ny;   % 网格点总数：21 × 21 = 441

grid_points = [X_grid(:), Y_grid(:)];

%% ============================================================
%  3. 将 9 个传感器映射到 21 × 21 网格索引
% ============================================================

sensor_idx = zeros(size(sensor_points,1),1);

for i = 1:size(sensor_points,1)

    x_coord = sensor_points(i,1);
    y_coord = sensor_points(i,2);

    col = x_coord + 1;
    row = y_coord + 1;

    sensor_idx(i) = sub2ind([ny, nx], row, col);

end

%% ============================================================
%  4. 根据 9 个传感器值估计接触中心、半径和幅值参数
% ============================================================

[cx_current, cy_current, r_current, amp_current, F_base_current, peak_current] = ...
    estimate_contact_params(current_sensor_values, sensor_points);

fprintf('\n当前接触参数估计结果：\n');
fprintf('接触中心：cx = %.4f mm, cy = %.4f mm\n', cx_current, cy_current);
fprintf('接触半径：r = %.4f mm\n', r_current);
fprintf('峰值力：peak = %.4f N\n', peak_current);
fprintf('基础力：F_base = %.4f N\n\n', F_base_current);

%% ============================================================
%  5. 基于加权质心和高斯模型生成 21 × 21 力场
% ============================================================

F_map_gaussian_current = F_base_current + amp_current * exp( ...
    -((X_grid - cx_current).^2 + (Y_grid - cy_current).^2) / ...
    (2 * r_current^2));

GaussianResultTable = table( ...
    grid_points(:,1), ...
    grid_points(:,2), ...
    F_map_gaussian_current(:), ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N'} ...
);

disp('加权质心 + 高斯拟合得到的当前力场结果：');
disp(GaussianResultTable);

%% ============================================================
%  6. 绘制高斯拟合力场
% ============================================================

figure;

surf(X_grid, Y_grid, F_map_gaussian_current, 'EdgeColor', 'none');
colormap jet;
colorbar;
ylabel(colorbar, 'Predicted force / N');

hold on;

scatter3( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    current_sensor_values, ...
    120, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

xlabel('x / mm');
ylabel('y / mm');
zlabel('Predicted force / N');

title('Force field by weighted centroid and Gaussian fitting');

axis tight;
grid on;
view(45, 30);

legend('Gaussian fitted force field', 'Sensor measurements');

%% ============================================================
%  7. 构造概率反演候选数据库
%
%  说明：
%  原 MLP 代码这里是构造训练数据集。
%  本代码不训练神经网络，而是构造候选压力场数据库。
%
%  每个候选样本包括：
%      X_data(k,:)：9 个传感器响应
%      Y_data(k,:)：对应的 21×21 力场
%
%  后续根据当前 9 个传感器值，
%  计算每个候选样本的后验概率。
% ============================================================

N = 3000;     % 候选样本数量，可根据需要增大，例如 10000

X_data = zeros(N, 9);
Y_data = zeros(N, n_grid);

candidate_cx = zeros(N,1);
candidate_cy = zeros(N,1);
candidate_sigma = zeros(N,1);
candidate_amp = zeros(N,1);
candidate_base = zeros(N,1);

for k = 1:N

    %% --------------------------------------------------------
    %  7.1 随机生成一个候选接触压力场
    %% --------------------------------------------------------

    cx_true = 20 * rand();              % 接触中心 x 坐标
    cy_true = 20 * rand();              % 接触中心 y 坐标

    sigma_true = 1.5 + 5.0 * rand();    % 接触区域大小
    amp_true = 800 + 6000 * rand();     % 压力峰值幅值
    F_base_true = 500 + 800 * rand();   % 基础力/背景力

    F_map_true = F_base_true + amp_true * exp( ...
        -((X_grid - cx_true).^2 + (Y_grid - cy_true).^2) / ...
        (2 * sigma_true^2));

    %% --------------------------------------------------------
    %  7.2 提取 9 个传感器位置处的响应
    %% --------------------------------------------------------

    sensor_values = F_map_true(sensor_idx);

    % 可选：加入传感器噪声
    % noise_level = 0.01;
    % sensor_values = sensor_values .* (1 + noise_level * randn(size(sensor_values)));

    %% --------------------------------------------------------
    %  7.3 根据 9 个传感器值重新估计高斯参数
    %
    %  这里与原 MLP 代码保持一致：
    %  不直接使用真实高斯场作为标签，
    %  而是使用 9 点估计后的高斯重建场作为候选标签。
    %% --------------------------------------------------------

    [cx_fit, cy_fit, r_fit, amp_fit, F_base_fit, ~] = ...
        estimate_contact_params(sensor_values, sensor_points);

    %% --------------------------------------------------------
    %  7.4 生成候选 21 × 21 力场
    %% --------------------------------------------------------

    F_map_fit = F_base_fit + amp_fit * exp( ...
        -((X_grid - cx_fit).^2 + (Y_grid - cy_fit).^2) / ...
        (2 * r_fit^2));

    X_data(k,:) = sensor_values(:)';
    Y_data(k,:) = F_map_fit(:)';

    candidate_cx(k) = cx_fit;
    candidate_cy(k) = cy_fit;
    candidate_sigma(k) = r_fit;
    candidate_amp(k) = amp_fit;
    candidate_base(k) = F_base_fit;

end

%% ============================================================
%  8. 设置概率反演噪声模型
%
%  说明：
%  概率反演不需要归一化训练。
%  它需要定义观测误差，即传感器测量值允许有多大偏差。
% ============================================================

relative_noise = 0.05;       % 5% 相对噪声，可调整为 0.01 ~ 0.10

sensor_scale = std(X_data, 0, 1);
fallback_scale = max(abs(X_data), [], 1);

sensor_scale(sensor_scale < eps) = fallback_scale(sensor_scale < eps);
sensor_scale(sensor_scale < eps) = 1;

sigma_sensor = relative_noise .* sensor_scale;

sigma_floor = max(1e-6, 0.001 .* max(abs(X_data), [], 1));
sigma_sensor = max(sigma_sensor, sigma_floor);

fprintf('\n概率反演噪声标准差：\n');
disp(sigma_sensor);

%% ============================================================
%  9. 设置候选样本先验概率
%
%  当前使用均匀先验：
%  每个候选压力场在观测前被认为同等可能。
% ============================================================

prior_probability = ones(N,1) / N;

%% ============================================================
%  10. 执行概率反演
%
%  输入：
%      当前 9 个传感器值
%
%  输出：
%      每个候选样本的后验概率
%      MAP 最优压力场
%      后验均值压力场
%      后验标准差压力场
% ============================================================

X_test = current_sensor_values(:)';

BayesResult = bayesian_inverse_predict( ...
    X_test, ...
    X_data, ...
    Y_data, ...
    sigma_sensor, ...
    prior_probability);

%% ============================================================
%  11. 输出概率反演参数结果
% ============================================================

map_idx = BayesResult.map_index;

fprintf('\n概率反演 MAP 最优结果：\n');
fprintf('候选编号：%d / %d\n', map_idx, N);
fprintf('接触中心：cx = %.4f mm, cy = %.4f mm\n', ...
    candidate_cx(map_idx), candidate_cy(map_idx));
fprintf('接触半径：r = %.4f mm\n', candidate_sigma(map_idx));
fprintf('峰值幅值：amp = %.4f N\n', candidate_amp(map_idx));
fprintf('基础力：F_base = %.4f N\n', candidate_base(map_idx));
fprintf('MAP 后验概率：%.8f\n', BayesResult.posterior_probability(map_idx));
fprintf('后验熵：%.6f\n', BayesResult.posterior_entropy);
fprintf('有效候选数：%.4f\n\n', BayesResult.effective_sample_size);

%% ============================================================
%  12. 输出后验概率最高的前 10 个候选样本
% ============================================================

top_k = min(10, N);

[sorted_prob, sorted_idx] = sort( ...
    BayesResult.posterior_probability, ...
    'descend');

TopTable = table( ...
    sorted_idx(1:top_k), ...
    candidate_cx(sorted_idx(1:top_k)), ...
    candidate_cy(sorted_idx(1:top_k)), ...
    candidate_sigma(sorted_idx(1:top_k)), ...
    candidate_amp(sorted_idx(1:top_k)), ...
    candidate_base(sorted_idx(1:top_k)), ...
    sorted_prob(1:top_k), ...
    'VariableNames', { ...
        'CandidateIndex', ...
        'cx_mm', ...
        'cy_mm', ...
        'r_mm', ...
        'amp_N', ...
        'F_base_N', ...
        'PosteriorProbability'});

disp('后验概率最高的候选样本：');
disp(TopTable);

%% ============================================================
%  13. 得到概率反演重建力场
%
%  MAP_Result：
%      后验概率最大的单个候选力场
%
%  Mean_Result：
%      所有候选力场按后验概率加权平均
%
%  Std_Result：
%      每个网格点的不确定性
% ============================================================

F_pred_map = BayesResult.map_pressure;
F_pred_mean = BayesResult.posterior_mean_pressure;
F_pred_std = BayesResult.posterior_std_pressure;

% 防止负值
F_pred_map(F_pred_map < 0) = 0;
F_pred_mean(F_pred_mean < 0) = 0;

F_map_pred = reshape(F_pred_map, ny, nx);
F_map_mean = reshape(F_pred_mean, ny, nx);
F_map_std = reshape(F_pred_std, ny, nx);

%% ============================================================
%  14. 输出概率反演预测结果表
% ============================================================

ResultTable = table( ...
    grid_points(:,1), ...
    grid_points(:,2), ...
    F_pred_map(:), ...
    F_pred_mean(:), ...
    F_pred_std(:), ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'MAPPredictedForce_N', ...
        'PosteriorMeanForce_N', ...
        'PosteriorStdForce_N'} ...
);

disp('概率反演重建得到的当前力场结果：');
disp(ResultTable);

%% ============================================================
%  15. 检查传感器位置处的预测值与输入值是否一致
% ============================================================

sensor_pred_values = F_map_pred(sensor_idx);

CompareTable = table( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    current_sensor_values, ...
    sensor_pred_values, ...
    current_sensor_values - sensor_pred_values, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'InputSensorForce_N', ...
        'BayesPredictedForce_N', ...
        'Error_N'} ...
);

disp('传感器位置处的输入值与概率反演预测值对比：');
disp(CompareTable);

%% ============================================================
%  16. 绘制概率反演 MAP 三维重建力场
% ============================================================

figure;

surf(X_grid, Y_grid, F_map_pred, 'EdgeColor', 'none');
colormap jet;
colorbar;
ylabel(colorbar, 'Predicted force / N');

hold on;

scatter3( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    current_sensor_values, ...
    120, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

xlabel('x / mm');
ylabel('y / mm');
zlabel('Predicted force / N');

title('Bayesian MAP reconstructed tactile force field');

axis tight;
grid on;
view(45, 30);

legend('Bayesian MAP surface', 'Sensor measurements');

%% ============================================================
%  17. 绘制概率反演二维等高线力场
% ============================================================

figure;

contourf(X_grid, Y_grid, F_map_pred, 30, 'LineColor', 'none');
hold on;

scatter( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    80, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('Bayesian MAP reconstructed force field');

%% ============================================================
%  18. 绘制后验均值力场和不确定性
% ============================================================

figure;

subplot(1,2,1);

contourf(X_grid, Y_grid, F_map_mean, 30, 'LineColor', 'none');
hold on;
scatter(sensor_points(:,1), sensor_points(:,2), 80, current_sensor_values, ...
    'filled', 'MarkerEdgeColor', 'k');
colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('Posterior mean force field');

subplot(1,2,2);

contourf(X_grid, Y_grid, F_map_std, 30, 'LineColor', 'none');
colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('Posterior standard deviation');

%% ============================================================
%  19. 对比高斯拟合结果与概率反演预测结果
% ============================================================

figure;

subplot(1,2,1);
contourf(X_grid, Y_grid, F_map_gaussian_current, 30, 'LineColor', 'none');
hold on;
scatter(sensor_points(:,1), sensor_points(:,2), 80, current_sensor_values, ...
    'filled', 'MarkerEdgeColor', 'k');
colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('Weighted centroid + Gaussian fit');

subplot(1,2,2);
contourf(X_grid, Y_grid, F_map_pred, 30, 'LineColor', 'none');
hold on;
scatter(sensor_points(:,1), sensor_points(:,2), 80, current_sensor_values, ...
    'filled', 'MarkerEdgeColor', 'k');
colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('Bayesian MAP prediction');
