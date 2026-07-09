clear;
clc;
close all;

%% ============================================================
%  MLP 压力场预测程序
%
%  功能：
%  1. 加载已经训练好的 MLP 模型；
%  2. 输入一组新的 9 点传感器值；
%  3. 预测 21×21 压力场；
%  4. 输出结果表；
%  5. 绘制三维图和二维等高线图；
%  6. 可选：与直接高斯拟合结果比较。
%
%  运行前：
%  请先运行 train_mlp_model.m，
%  生成 trained_pressure_model.mat。
%% ============================================================

%% 1. 加载模型

model_file = 'trained_pressure_model.mat';

if ~isfile(model_file)
    error(['未找到模型文件：', model_file, ...
        newline, ...
        '请先运行 train_mlp_model.m 完成训练。']);
end

load(model_file);

fprintf('已加载模型：%s\n', model_file);
fprintf('模型验证 RMSE：%.4f N\n', val_RMSE);

%% 2. 输入当前 9 个传感器值

% 传感器顺序必须与训练时 sensor_points 的顺序完全一致：
%
% (0,0), (10,0), (20,0),
% (0,10), (10,10), (20,10),
% (0,20), (10,20), (20,20)

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

if numel(current_sensor_values) ~= 9
    error('current_sensor_values 必须包含 9 个传感器值。');
end

if any(~isfinite(current_sensor_values))
    error('传感器值中包含 NaN 或 Inf，请检查输入。');
end

%% 3. 归一化当前输入

X_test = current_sensor_values(:)';

X_test_norm = (X_test - mu_X) ./ std_X;

%% 4. 使用训练好的 MLP 预测

F_pred_norm = predict(net, X_test_norm);

F_pred = F_pred_norm .* std_Y + mu_Y;

% 删除不符合物理意义的负值
F_pred(F_pred < 0) = 0;

F_map_pred = reshape(F_pred, ny, nx);

%% 5. 输出预测结果表

ResultTable = table( ...
    grid_points(:, 1), ...
    grid_points(:, 2), ...
    F_pred(:), ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N'} ...
);

disp('MLP 预测得到的 21×21 力场：');
disp(ResultTable);

writetable(ResultTable, 'MLP_prediction_result.csv');

fprintf('预测结果已保存：MLP_prediction_result.csv\n');

%% 6. 检查传感器位置处的预测值

sensor_pred_values = F_map_pred(sensor_idx);

CompareTable = table( ...
    sensor_points(:, 1), ...
    sensor_points(:, 2), ...
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

disp('传感器输入值与 MLP 预测值对比：');
disp(CompareTable);

writetable(CompareTable, 'sensor_prediction_comparison.csv');

%% 7. 绘制 MLP 三维压力场

figure('Name', 'MLP 3D Force Field');

surf(X_grid, Y_grid, F_map_pred, ...
    'EdgeColor', 'none');

colormap jet;
cb = colorbar;
ylabel(cb, 'Predicted force / N');

hold on;

scatter3( ...
    sensor_points(:, 1), ...
    sensor_points(:, 2), ...
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

legend( ...
    'MLP reconstructed surface', ...
    'Sensor measurements', ...
    'Location', 'best');

%% 8. 绘制 MLP 二维等高线图

figure('Name', 'MLP 2D Force Field');

contourf( ...
    X_grid, ...
    Y_grid, ...
    F_map_pred, ...
    30, ...
    'LineColor', 'none');

hold on;

scatter( ...
    sensor_points(:, 1), ...
    sensor_points(:, 2), ...
    80, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

colorbar;

xlabel('x / mm');
ylabel('y / mm');

axis equal;
axis tight;

title('MLP reconstructed force field');

%% 9. 可选：计算直接高斯拟合结果，并与 MLP 比较

[cx_current, cy_current, r_current, ...
 amp_current, F_base_current, peak_current] = ...
    estimate_contact_params( ...
        current_sensor_values, ...
        sensor_points);

F_map_gaussian = F_base_current + ...
    amp_current * exp( ...
    -((X_grid - cx_current).^2 + ...
      (Y_grid - cy_current).^2) / ...
    (2 * r_current^2));

fprintf('\n直接高斯拟合参数：\n');
fprintf('接触中心：cx = %.4f mm, cy = %.4f mm\n', ...
    cx_current, cy_current);
fprintf('接触半径：r = %.4f mm\n', r_current);
fprintf('峰值：%.4f N\n', peak_current);
fprintf('基础值：%.4f N\n', F_base_current);

figure('Name', 'Gaussian versus MLP');

tiledlayout(1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile;

contourf( ...
    X_grid, ...
    Y_grid, ...
    F_map_gaussian, ...
    30, ...
    'LineColor', 'none');

hold on;

scatter( ...
    sensor_points(:, 1), ...
    sensor_points(:, 2), ...
    80, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
axis tight;
title('Weighted centroid + Gaussian fit');

nexttile;

contourf( ...
    X_grid, ...
    Y_grid, ...
    F_map_pred, ...
    30, ...
    'LineColor', 'none');

hold on;

scatter( ...
    sensor_points(:, 1), ...
    sensor_points(:, 2), ...
    80, ...
    current_sensor_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
axis tight;
title('MLP prediction');
