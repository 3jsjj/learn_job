clear;
clc;
close all;

rng(0);

%% ============================================================
%  MLP 训练程序
%
%  功能：
%  1. 使用高斯模型生成大量模拟训练样本；
%  2. 将 9 个传感器值作为网络输入；
%  3. 将 21×21 力场作为网络输出；
%  4. 训练并验证 MLP；
%  5. 保存网络和归一化参数。
%
%  注意：
%  当前训练标签仍由高斯模型生成。
%  后续若使用 Abaqus 数据，只需要替换“训练数据生成”部分。
%% ============================================================

%% 1. 定义传感器位置

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

%% 2. 构建 21×21 网格

x_int = 0:1:20;
y_int = 0:1:20;

[X_grid, Y_grid] = meshgrid(x_int, y_int);

nx = length(x_int);
ny = length(y_int);
n_grid = nx * ny;

grid_points = [X_grid(:), Y_grid(:)];

%% 3. 将 9 个传感器位置映射到网格索引

sensor_idx = zeros(size(sensor_points, 1), 1);

for i = 1:size(sensor_points, 1)

    x_coord = sensor_points(i, 1);
    y_coord = sensor_points(i, 2);

    col = x_coord + 1;
    row = y_coord + 1;

    sensor_idx(i) = sub2ind([ny, nx], row, col);

end

%% 4. 生成 MLP 训练数据

N = 3000;

X_data = zeros(N, 9);
Y_data = zeros(N, n_grid);

fprintf('开始生成训练数据，共 %d 组。\n', N);

for k = 1:N

    % ---------------------------------------------------------
    % 随机生成一组高斯接触参数
    % ---------------------------------------------------------

    cx_true = 20 * rand();
    cy_true = 20 * rand();

    sigma_true = 1.5 + 5.0 * rand();
    amp_true = 800 + 6000 * rand();
    F_base_true = 500 + 800 * rand();

    % ---------------------------------------------------------
    % 生成完整的 21×21 高斯压力场
    % ---------------------------------------------------------

    F_map_true = F_base_true + amp_true * exp( ...
        -((X_grid - cx_true).^2 + (Y_grid - cy_true).^2) / ...
        (2 * sigma_true^2));

    % ---------------------------------------------------------
    % 提取 9 个传感器位置的值，作为 MLP 输入
    % ---------------------------------------------------------

    sensor_values = F_map_true(sensor_idx);

    % 可选：加入传感器噪声
    % noise_level = 0.01;
    % sensor_values = sensor_values .* ...
    %     (1 + noise_level * randn(size(sensor_values)));

    % ---------------------------------------------------------
    % 根据 9 点数据重新估计高斯参数
    % ---------------------------------------------------------

    [cx_fit, cy_fit, r_fit, amp_fit, F_base_fit, ~] = ...
        estimate_contact_params(sensor_values, sensor_points);

    % ---------------------------------------------------------
    % 用估计参数生成训练标签
    % ---------------------------------------------------------

    F_map_fit = F_base_fit + amp_fit * exp( ...
        -((X_grid - cx_fit).^2 + (Y_grid - cy_fit).^2) / ...
        (2 * r_fit^2));

    X_data(k, :) = sensor_values(:)';
    Y_data(k, :) = F_map_fit(:)';

end

fprintf('训练数据生成完成。\n');

%% 5. 数据归一化

mu_X = mean(X_data, 1);
std_X = std(X_data, 0, 1) + eps;

mu_Y = mean(Y_data, 1);
std_Y = std(Y_data, 0, 1) + eps;

X_norm = (X_data - mu_X) ./ std_X;
Y_norm = (Y_data - mu_Y) ./ std_Y;

%% 6. 划分训练集和验证集

idx = randperm(N);

N_train = round(0.8 * N);

train_idx = idx(1:N_train);
val_idx = idx(N_train + 1:end);

X_train = X_norm(train_idx, :);
Y_train = Y_norm(train_idx, :);

X_val = X_norm(val_idx, :);
Y_val = Y_norm(val_idx, :);

fprintf('训练集：%d 组。\n', size(X_train, 1));
fprintf('验证集：%d 组。\n', size(X_val, 1));

%% 7. 定义 MLP 网络

layers = [
    featureInputLayer(9, ...
        'Normalization', 'none', ...
        'Name', 'input')

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

%% 8. 设置训练参数

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

%% 9. 训练 MLP

fprintf('开始训练 MLP。\n');

net = trainNetwork(X_train, Y_train, layers, options);

fprintf('MLP 训练完成。\n');

%% 10. 验证模型

Y_val_pred_norm = predict(net, X_val);

Y_val_pred = Y_val_pred_norm .* std_Y + mu_Y;
Y_val_true = Y_val .* std_Y + mu_Y;

val_RMSE = sqrt(mean((Y_val_pred(:) - Y_val_true(:)).^2));
val_MAE = mean(abs(Y_val_pred(:) - Y_val_true(:)));

fprintf('\n模型验证结果：\n');
fprintf('Validation RMSE = %.4f N\n', val_RMSE);
fprintf('Validation MAE  = %.4f N\n', val_MAE);

%% 11. 保存训练模型

model_file = 'trained_pressure_model.mat';

save(model_file, ...
    'net', ...
    'mu_X', ...
    'std_X', ...
    'mu_Y', ...
    'std_Y', ...
    'sensor_points', ...
    'sensor_idx', ...
    'x_int', ...
    'y_int', ...
    'X_grid', ...
    'Y_grid', ...
    'grid_points', ...
    'nx', ...
    'ny', ...
    'n_grid', ...
    'val_RMSE', ...
    'val_MAE');

fprintf('\n模型已保存：%s\n', model_file);
fprintf('以后运行 predict_force_field.m 即可直接预测，无须重新训练。\n');
