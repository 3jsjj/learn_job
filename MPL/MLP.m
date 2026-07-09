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
%
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
%
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
%  7. 构造 MLP 训练数据集
%
%
%  1. 定义传感器位置和当前测量值
%  2. 构建 21 × 21 插值网格
%  3. 将 9 个传感器映射到 21 × 21 网格索引
%
% ============================================================

N = 3000;     % 训练样本数量，可根据需要增大，例如 10000

X_data = zeros(N, 9);
Y_data = zeros(N, n_grid);

for k = 1:N

    %% --------------------------------------------------------
%  7. 构造 MLP 训练数据集
    %
    %% --------------------------------------------------------

    cx_true = 20 * rand();              % 接触中心 x 坐标
    cy_true = 20 * rand();              % 接触中心 y 坐标

    sigma_true = 1.5 + 5.0 * rand();    % 高斯分布标准差，控制接触区域大小
    amp_true = 800 + 6000 * rand();     % 压力峰值幅值
    F_base_true = 500 + 800 * rand();   % 基础力/背景力

    F_map_true = F_base_true + amp_true * exp( ...
        -((X_grid - cx_true).^2 + (Y_grid - cy_true).^2) / ...
        (2 * sigma_true^2));

    %% --------------------------------------------------------
%  7. 构造 MLP 训练数据集
    %% --------------------------------------------------------

    sensor_values = F_map_true(sensor_idx);

    % noise_level = 0.01;
    % sensor_values = sensor_values .* (1 + noise_level * randn(size(sensor_values)));

    %% --------------------------------------------------------
%  7. 构造 MLP 训练数据集
    %
    %% --------------------------------------------------------

    [cx_fit, cy_fit, r_fit, amp_fit, F_base_fit, ~] = ...
        estimate_contact_params(sensor_values, sensor_points);

    %% --------------------------------------------------------
%  7. 构造 MLP 训练数据集
    %% --------------------------------------------------------

    F_map_fit = F_base_fit + amp_fit * exp( ...
        -((X_grid - cx_fit).^2 + (Y_grid - cy_fit).^2) / ...
        (2 * r_fit^2));

    X_data(k,:) = sensor_values(:)';
    Y_data(k,:) = F_map_fit(:)';

end

%% ============================================================
%  8. 数据归一化
% ============================================================

mu_X = mean(X_data, 1);
std_X = std(X_data, 0, 1) + eps;

mu_Y = mean(Y_data, 1);
std_Y = std(Y_data, 0, 1) + eps;

X_norm = (X_data - mu_X) ./ std_X;
Y_norm = (Y_data - mu_Y) ./ std_Y;

%% ============================================================
%  9. 划分训练集和验证集
% ============================================================

idx = randperm(N);

N_train = round(0.8 * N);

train_idx = idx(1:N_train);
val_idx = idx(N_train+1:end);

X_train = X_norm(train_idx,:);
Y_train = Y_norm(train_idx,:);

X_val = X_norm(val_idx,:);
Y_val = Y_norm(val_idx,:);

%% ============================================================
%  10. 定义 MLP 网络结构
%
% ============================================================

layers = [
    featureInputLayer(9, 'Normalization', 'none', 'Name', 'input')

    fullyConnectedLayer(128, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    dropoutLayer(0.1, 'Name', 'dropout1')

    fullyConnectedLayer(256, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    dropoutLayer(0.1, 'Name', 'dropout2')

    fullyConnectedLayer(512, 'Name', 'fc3')
    reluLayer('Name', 'relu3')

    fullyConnectedLayer(441, 'Name', 'output')

    regressionLayer('Name', 'regression')
];

%% ============================================================
%  11. 设置训练参数
% ============================================================

options = trainingOptions('adam', ...
    'MaxEpochs', 80, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 1e-3, ...
    'L2Regularization', 1e-4, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, Y_val}, ...
    'ValidationFrequency', 30, ...
    'ValidationPatience', 10, ...
    'Plots', 'training-progress', ...
    'Verbose', false);

%% ============================================================
%  12. 训练 MLP 网络
% ============================================================

net = trainNetwork(X_train, Y_train, layers, options);

%% ============================================================
%  13. 验证集预测与误差评估
% ============================================================

Y_val_pred_norm = predict(net, X_val);

Y_val_pred = Y_val_pred_norm .* std_Y + mu_Y;
Y_val_true = Y_val .* std_Y + mu_Y;

val_RMSE = sqrt(mean((Y_val_pred(:) - Y_val_true(:)).^2));

fprintf('Validation RMSE = %.4f N\n', val_RMSE);

%% ============================================================
%  14. 使用当前 9 个传感器值重建 21 × 21 力场
% ============================================================

X_test = current_sensor_values(:)';

X_test_norm = (X_test - mu_X) ./ std_X;

F_pred_norm = predict(net, X_test_norm);

F_pred = F_pred_norm .* std_Y + mu_Y;

F_pred(F_pred < 0) = 0;

F_map_pred = reshape(F_pred, ny, nx);

%% ============================================================
%  15. 输出 MLP 预测结果表
% ============================================================

ResultTable = table( ...
    grid_points(:,1), ...
    grid_points(:,2), ...
    F_pred(:), ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N'} ...
);

disp('MLP 重建得到的当前力场结果：');
disp(ResultTable);

%% ============================================================
%  16. 检查传感器位置处的预测值与输入值是否一致
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
        'MLPPredictedForce_N', ...
        'Error_N'} ...
);

disp('传感器位置处的输入值与 MLP 预测值对比：');
disp(CompareTable);

%% ============================================================
%  17. 绘制 MLP 三维重建力场
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

title('MLP reconstructed tactile force field');

axis tight;
grid on;
view(45, 30);

legend('MLP reconstructed surface', 'Sensor measurements');

%% ============================================================
%  18. 绘制 MLP 二维等高线力场
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
title('MLP reconstructed force field');

%% ============================================================
%  19. 对比高斯拟合结果与 MLP 预测结果
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
title('MLP prediction');
