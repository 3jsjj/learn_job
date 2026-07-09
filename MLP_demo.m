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
n_grid = nx * ny;   % 441

grid_points = [X_grid(:), Y_grid(:)];

%% ============================================================
%  3. 找到 9 个传感器在 21×21 网格中的索引
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
%  4. 生成 MLP 训练数据
%
%  输入 X_data: N × 9
%  输出 Y_data: N × 441
%
%  这里用随机高斯接触力场做 demo。
%  正式使用时，可以把这里换成 Abaqus / 实验数据。
% ============================================================

N = 3000;     % 训练样本数量，可以改大，比如 10000

X_data = zeros(N, 9);
Y_data = zeros(N, n_grid);

for k = 1:N

    F_map = zeros(ny, nx);

    % 每个样本随机生成 1~3 个接触区域
    n_contact = randi([1, 3]);

    for c = 1:n_contact

        cx = 20 * rand();              % 接触中心 x
        cy = 20 * rand();              % 接触中心 y
        amp = 800 + 6000 * rand();     % 接触峰值力
        sigma = 1.5 + 3.5 * rand();    % 接触扩散半径

        F_map = F_map + amp * exp( ...
            -((X_grid - cx).^2 + (Y_grid - cy).^2) / (2 * sigma^2));

    end

    % 9 个传感器位置的力值作为输入
    sensor_values = F_map(sensor_idx);

    % 可选：加入少量噪声，模拟测量误差
    noise_level = 0.01;
    sensor_values = sensor_values .* (1 + noise_level * randn(size(sensor_values)));

    X_data(k,:) = sensor_values(:)';
    Y_data(k,:) = F_map(:)';

end

%% ============================================================
%  5. 数据归一化
% ============================================================

mu_X = mean(X_data, 1);
std_X = std(X_data, 0, 1) + eps;

mu_Y = mean(Y_data, 1);
std_Y = std(Y_data, 0, 1) + eps;

X_norm = (X_data - mu_X) ./ std_X;
Y_norm = (Y_data - mu_Y) ./ std_Y;

%% ============================================================
%  6. 划分训练集和验证集
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
%  7. 建立 MLP 网络
%
%  输入：9
%  输出：441
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
%  8. 训练参数
% ============================================================

options = trainingOptions('adam', ...
    'MaxEpochs', 80, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, Y_val}, ...
    'ValidationFrequency', 30, ...
    'Plots', 'training-progress', ...
    'Verbose', false);

%% ============================================================
%  9. 训练 MLP
% ============================================================

net = trainNetwork(X_train, Y_train, layers, options);

%% ============================================================
%  10. 验证集误差
% ============================================================

Y_val_pred_norm = predict(net, X_val);

Y_val_pred = Y_val_pred_norm .* std_Y + mu_Y;
Y_val_true = Y_val .* std_Y + mu_Y;

val_RMSE = sqrt(mean((Y_val_pred(:) - Y_val_true(:)).^2));

fprintf('Validation RMSE = %.4f N\n', val_RMSE);

%% ============================================================
%  11. 输入当前 9 个霍尔力值，预测 21×21 力场
% ============================================================

X_test = current_sensor_values(:)';

X_test_norm = (X_test - mu_X) ./ std_X;

F_pred_norm = predict(net, X_test_norm);

F_pred = F_pred_norm .* std_Y + mu_Y;

F_map_pred = reshape(F_pred, ny, nx);

%% ============================================================
%  12. 输出所有整数坐标点的力数据
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
%  13. 检查传感器位置预测值和输入值是否一致
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
%  14. 绘制 MLP 重建力场
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
%  15. 二维等高线图
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