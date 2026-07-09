clc;
clear;
close all;

%% 1. 真实霍尔测点坐标 (mm)
X_meas = [
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

%% 2. 霍尔芯片计算出的力 (N)
F_meas = [
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

%% 3. 创建连续预测网格，用于画图
x = linspace(0, 20, 50);
y = linspace(0, 20, 50);

[xx, yy] = meshgrid(x, y);
X_virtual = [xx(:), yy(:)];

%% 4. 建立 Gaussian Process Regression 模型
gprMdl = fitrgp( ...
    X_meas, ...
    F_meas, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...
    'Standardize', true, ...
    'Sigma', 0.05, ...
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact');

%% 5. 预测连续网格力分布
[F_pred, F_std] = predict(gprMdl, X_virtual);

F_map = reshape(F_pred, size(xx));
Std_map = reshape(F_std, size(xx));

%% 6. 绘制力分布图
figure;

contourf(xx, yy, F_map, 30, 'LineColor', 'none');
hold on;

scatter( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    80, ...
    F_meas, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

colorbar;
xlabel('x / mm');
ylabel('y / mm');
axis equal;
title('GP reconstructed force field');

%% 7. 输出所有整数坐标点的力数据

x_int = 0:1:20;
y_int = 0:1:20;

[xx_int, yy_int] = meshgrid(x_int, y_int);

X_int = [xx_int(:), yy_int(:)];

[F_int, F_int_std] = predict(gprMdl, X_int);

ResultTable = table( ...
    X_int(:,1), ...
    X_int(:,2), ...
    F_int, ...
    F_int_std, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N', ...
        'Uncertainty_N'} ...
);

disp('整数坐标点力数据：');
disp(ResultTable);

%% 8. 保存整数坐标点力数据到 Excel
writetable(ResultTable, 'integer_coordinate_force_result.xlsx');

disp('整数坐标点力数据已保存为 integer_coordinate_force_result.xlsx');

%% 9. 输出学习后的模型参数
disp('Learned kernel parameters:');
disp(gprMdl.KernelInformation);

disp('Sigma noise:');
disp(gprMdl.Sigma);