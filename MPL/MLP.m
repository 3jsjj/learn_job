clear;
clc;
close all;

rng(0);

%% ============================================================
%  1. 传感器位置和当前测量力值
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
%  2. 建立 21 × 21 整数坐标网格
% ============================================================

x_int = 0:1:20;
y_int = 0:1:20;

[X_grid, Y_grid] = meshgrid(x_int, y_int);

nx = length(x_int);
ny = length(y_int);
n_grid = nx * ny;   % 441 个整数坐标点

grid_points = [X_grid(:), Y_grid(:)];

%% ============================================================
%  3. 找到 9 个传感器在 21 × 21 网格中的索引
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
%  4. 根据当前 9 个传感器力值估计接触参数
%
%  接触中心：加权质心法
%  接触半径：不手动指定，通过误差最小化自动拟合
%  峰值力：在拟合半径下，通过高斯最小二乘拟合
% ============================================================

[cx_current, cy_current, r_current, amp_current, F_base_current, peak_current] = ...
    estimate_contact_params(current_sensor_values, sensor_points);

fprintf('\n当前输入的接触参数估计结果：\n');
fprintf('接触中心 cx = %.4f mm, cy = %.4f mm\n', cx_current, cy_current);
fprintf('自动拟合接触半径 r = %.4f mm\n', r_current);
fprintf('高斯拟合峰值力 peak = %.4f N\n', peak_current);
fprintf('基础力 F_base = %.4f N\n\n', F_base_current);

%% ============================================================
%  5. 根据估计参数生成当前 21 × 21 高斯拟合力场
%
%  这一步不是 MLP，而是用：
%  加权质心 + 自动半径拟合 + 高斯峰值拟合
%  直接得到一个参数化力场，用于和 MLP 结果对比
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

disp('基于加权质心 + 高斯拟合的整数坐标点力数据：');
disp(GaussianResultTable);

%% ============================================================
%  6. 绘制当前高斯拟合力场
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
%  7. 生成 MLP 训练数据
%
%  输入 X_data: N × 9
%  输出 Y_data: N × 441
%
%  说明：
%  这里仍然需要生成大量训练样本。
%  但是每个样本的输出力场不再直接使用随机参数，
%  而是先由 9 个传感器值估计：
%      1. 接触中心：加权质心
%      2. 接触半径：自动拟合
%      3. 峰值力：高斯拟合
%  然后生成对应的 21 × 21 力场作为训练标签。
%
%  正式使用时，建议把这里替换成 Abaqus 或实验数据。
% ============================================================

N = 3000;     % 训练样本数量，可改大，例如 10000

X_data = zeros(N, 9);
Y_data = zeros(N, n_grid);

for k = 1:N

    %% --------------------------------------------------------
    %  7.1 先生成一个隐藏的单峰接触样本
    %
    %  注意：
    %  这里的随机参数只是为了构造 demo 训练样本。
    %  最终训练标签不是直接用这些随机参数，
    %  而是通过 9 个传感器值重新估计接触中心、半径和峰值力。
    %% --------------------------------------------------------

    cx_true = 20 * rand();              % 隐藏真实接触中心 x
    cy_true = 20 * rand();              % 隐藏真实接触中心 y

    sigma_true = 1.5 + 5.0 * rand();    % 隐藏真实扩散半径，仅用于生成样本
    amp_true = 800 + 6000 * rand();     % 隐藏真实峰值幅值
    F_base_true = 500 + 800 * rand();   % 隐藏基础力

    F_map_true = F_base_true + amp_true * exp( ...
        -((X_grid - cx_true).^2 + (Y_grid - cy_true).^2) / ...
        (2 * sigma_true^2));

    %% --------------------------------------------------------
    %  7.2 提取 9 个传感器位置的力值作为 MLP 输入
    %% --------------------------------------------------------

    sensor_values = F_map_true(sensor_idx);

    % 可选：加入少量噪声，模拟测量误差
    % noise_level = 0.01;
    % sensor_values = sensor_values .* (1 + noise_level * randn(size(sensor_values)));

    %% --------------------------------------------------------
    %  7.3 根据 9 个传感器值估计接触参数
    %
    %  接触中心：加权质心
    %  接触半径：自动拟合
    %  峰值力：高斯拟合
    %% --------------------------------------------------------

    [cx_fit, cy_fit, r_fit, amp_fit, F_base_fit, ~] = ...
        estimate_contact_params(sensor_values, sensor_points);

    %% --------------------------------------------------------
    %  7.4 用估计参数生成 21 × 21 训练标签
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
%  10. 建立 MLP 神经网络
%
%  输入：9 个传感器力值
%  输出：441 个整数坐标点力值
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
%  13. 计算验证集误差
% ============================================================

Y_val_pred_norm = predict(net, X_val);

Y_val_pred = Y_val_pred_norm .* std_Y + mu_Y;
Y_val_true = Y_val .* std_Y + mu_Y;

val_RMSE = sqrt(mean((Y_val_pred(:) - Y_val_true(:)).^2));

fprintf('Validation RMSE = %.4f N\n', val_RMSE);

%% ============================================================
%  14. 输入当前 9 个霍尔传感器力值，预测 21 × 21 力场
% ============================================================

X_test = current_sensor_values(:)';

X_test_norm = (X_test - mu_X) ./ std_X;

F_pred_norm = predict(net, X_test_norm);

F_pred = F_pred_norm .* std_Y + mu_Y;

% 防止出现负力
F_pred(F_pred < 0) = 0;

F_map_pred = reshape(F_pred, ny, nx);

%% ============================================================
%  15. 输出所有整数坐标点的 MLP 预测力数据
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

disp('MLP预测的整数坐标点力数据：');
disp(ResultTable);

%% ============================================================
%  16. 检查传感器位置预测值和输入值是否一致
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

disp('传感器位置拟合对比：');
disp(CompareTable);

%% ============================================================
%  17. 绘制 MLP 重建的三维力场
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
%  18. 绘制 MLP 二维等高线力场图
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
%  19. 对比：高斯拟合结果和 MLP 结果
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