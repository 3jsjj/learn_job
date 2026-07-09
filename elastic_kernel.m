clear;
clc;
close all;

%% 1. 真实测点坐标 (mm)
real_points = [
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

%% 2. 真实测点力值 (N)
real_values = [
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

%% 3. 线弹性材料参数
E = 100000.0;     % Pa
nu = 0.49;

G = E / (2 * (1 + nu));

lambda_elastic = 8.0;   % 扩散长度，单位与坐标一致
eps_val = 1e-6;

%% 4. 创建连续虚拟网格，用于画图
x = linspace(0, 20, 50);
y = linspace(0, 20, 50);

[X, Y] = meshgrid(x, y);

virtual_points = [X(:), Y(:)];

%% 5. 计算连续网格上的弹性扩散力分布
virtual_values = zeros(size(virtual_points, 1), 1);

for i = 1:size(virtual_points, 1)

    vp = virtual_points(i, :);

    distances = sqrt( ...
        (real_points(:,1) - vp(1)).^2 + ...
        (real_points(:,2) - vp(2)).^2 );

    weights = exp(-distances / lambda_elastic) ./ ...
              (G * (distances + eps_val));

    weights = weights / sum(weights);

    virtual_values(i) = sum(weights .* real_values);

end

Z = reshape(virtual_values, size(X));

%% 6. 三维可视化
figure;

surf(X, Y, Z, 'EdgeColor', 'none');
colormap jet;
colorbar;
ylabel(colorbar, 'Predicted force / N');

hold on;

scatter3( ...
    real_points(:,1), ...
    real_points(:,2), ...
    real_values, ...
    120, ...
    real_values, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

xlabel('x / mm');
ylabel('y / mm');
zlabel('Predicted force / N');

title('3D virtual tactile field by elastic diffusion');

axis tight;
grid on;
view(45, 30);

legend('Interpolated surface', 'Real sensors');

%% 7. 输出所有整数坐标点的力数据

x_int = 0:1:20;
y_int = 0:1:20;

[X_int_grid, Y_int_grid] = meshgrid(x_int, y_int);

integer_points = [X_int_grid(:), Y_int_grid(:)];

integer_values = zeros(size(integer_points, 1), 1);

for i = 1:size(integer_points, 1)

    ip = integer_points(i, :);

    distances = sqrt( ...
        (real_points(:,1) - ip(1)).^2 + ...
        (real_points(:,2) - ip(2)).^2 );

    weights = exp(-distances / lambda_elastic) ./ ...
              (G * (distances + eps_val));

    weights = weights / sum(weights);

    integer_values(i) = sum(weights .* real_values);

end

ResultTable = table( ...
    integer_points(:,1), ...
    integer_points(:,2), ...
    integer_values, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N'} ...
);

disp('整数坐标点力数据：');
disp(ResultTable);