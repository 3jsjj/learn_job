clear;
clc;
close all;

%% 1. 真实测点坐标 (mm)
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

%% 2. 真实测点力值 (N)
b = [
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

%% 3. 建立整数网格 21 × 21
x_int = 0:1:20;
y_int = 0:1:20;

[X_grid, Y_grid] = meshgrid(x_int, y_int);

grid_points = [X_grid(:), Y_grid(:)];

n_grid = size(grid_points, 1);     % 441
n_sensor = size(sensor_points, 1); % 9

%% 4. 构造测量矩阵 A
% A 表示每个网格点对每个传感器测量值的贡献
% 这里使用 RBF 形式的空间响应函数

sigma = 4.0;  % 空间影响半径，单位 mm

A = zeros(n_sensor, n_grid);

for i = 1:n_sensor
    for j = 1:n_grid

        d = norm(sensor_points(i,:) - grid_points(j,:));

        A(i,j) = exp(-(d^2) / (2 * sigma^2));

    end
end

%% 5. 构造二维 Laplacian 正则化矩阵 L
% L 用来约束重建结果平滑

nx = length(x_int);
ny = length(y_int);

N = nx * ny;

L = sparse(N, N);

for row = 1:ny
    for col = 1:nx

        idx = sub2ind([ny, nx], row, col);

        L(idx, idx) = -4;

        if row > 1
            idx_up = sub2ind([ny, nx], row - 1, col);
            L(idx, idx_up) = 1;
        end

        if row < ny
            idx_down = sub2ind([ny, nx], row + 1, col);
            L(idx, idx_down) = 1;
        end

        if col > 1
            idx_left = sub2ind([ny, nx], row, col - 1);
            L(idx, idx_left) = 1;
        end

        if col < nx
            idx_right = sub2ind([ny, nx], row, col + 1);
            L(idx, idx_right) = 1;
        end
    end
end

%% 6. Tikhonov 正则化求解
lambda = 0.1;

f = (A' * A + lambda * (L' * L)) \ (A' * b);

F_map = reshape(f, ny, nx);

%% 7. 绘制重建力场
figure;

surf(X_grid, Y_grid, F_map, 'EdgeColor', 'none');
colormap jet;
colorbar;
ylabel(colorbar, 'Predicted force / N');

hold on;

scatter3( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    b, ...
    120, ...
    b, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

xlabel('x / mm');
ylabel('y / mm');
zlabel('Predicted force / N');

title('Tikhonov reconstructed tactile field');

axis tight;
grid on;
view(45, 30);

legend('Tikhonov surface', 'Real sensors');

%% 8. 输出所有整数坐标点的力数据

ResultTable = table( ...
    grid_points(:,1), ...
    grid_points(:,2), ...
    f, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N'} ...
);

disp('整数坐标点力数据：');
disp(ResultTable);

%% 9. 检查在传感器位置的拟合效果
b_fit = A * f;

CompareTable = table( ...
    sensor_points(:,1), ...
    sensor_points(:,2), ...
    b, ...
    b_fit, ...
    b - b_fit, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'MeasuredForce_N', ...
        'FittedForce_N', ...
        'Error_N'} ...
);

disp('传感器位置拟合效果：');
disp(CompareTable);