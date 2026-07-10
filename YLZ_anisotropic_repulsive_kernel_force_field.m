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

%% 3. YLZ 各向异性排斥核参数
% 径向部分按以下有限作用范围形式构造：
%
%   u_R(r) = 2*epsilon0,                              0 <= r <= r_min
%
%   u_R(r) = 2*epsilon0*cos( pi*(r-r_min)
%                    / (2*(r_c-r_min)) )^(2*zeta),   r_min < r < r_c
%
%   u_R(r) = 0,                                      r >= r_c
%
% 方向部分：
%
%   phi = 1 + mu*(a - 1)
%
%   a = (n_i x r_hat)*(n_j x r_hat)
%       + beta*(n_i-n_j)·r_hat - beta^2
%
% 最终核：
%
%   K_ij = u_R(r_ij)*phi_ij
%
% 为了将其作为非负插值权重，程序会执行 phi = max(phi, 0)。

params.epsilon0 = 1.0;     % 势函数幅值；归一化插值时公共幅值会被约掉
params.r_min    = 2.0;     % 近场半径 (mm)
params.r_c      = 12.0;    % 截止半径 (mm)，必须大于 r_min
params.zeta     = 2.0;     % 径向衰减陡峭程度
params.mu       = 0.8;     % 各向异性强度
params.beta     = 0.2;     % 方向偏置参数

params.r_eps        = 1e-12;  % 避免 r = 0 时计算单位方向向量出现除零
params.weight_tol   = 1e-12;  % 判断权重和是否为零
params.coincide_tol = 1e-10;  % 判断预测点是否与测点重合

if params.r_c <= params.r_min
    error('参数错误：r_c 必须大于 r_min。');
end

%% 4. 设置真实测点和虚拟点的方向向量
% 每个方向向量均定义在 x-y 平面内，并应当为单位向量。
%
% 当前示例令所有测点和虚拟点方向均沿 x 轴：
%
%   n = [1, 0]
%
% 这时 beta*(n_i-n_j)·r_hat 项为 0，但交叉乘积项仍然会产生
% 与传播方向有关的各向异性。
%
% 若每个霍尔测点具有不同方向，可分别修改 theta_meas_deg。
% 若虚拟位置具有空间变化的材料方向，也可单独构造 N_virtual。

theta_meas_deg = zeros(size(X_meas, 1), 1);

N_meas = [
    cosd(theta_meas_deg), ...
    sind(theta_meas_deg)
];

%% 5. 创建连续预测网格，用于画图
x = linspace(0, 20, 50);
y = linspace(0, 20, 50);

[xx, yy] = meshgrid(x, y);
X_virtual = [xx(:), yy(:)];

% 当前假设所有虚拟点的方向也沿 x 轴。
theta_virtual_deg = zeros(size(X_virtual, 1), 1);

N_virtual = [
    cosd(theta_virtual_deg), ...
    sind(theta_virtual_deg)
];

%% 6. 使用 YLZ 各向异性排斥核预测连续力场
% 多测点重构采用归一化核加权：
%
%          sum_j K_ij * F_j
%   F_i = -------------------
%              sum_j K_ij
%
% 这是一种方向相关的核插值，不是 Gaussian Process。
% 因此它不会直接给出严格的概率预测标准差。
%
% F_dispersion 是邻近测点力值的加权局部离散度，
% 只能用于表示局部测点之间的一致程度，不能视为 GP 不确定度。

[F_pred, F_dispersion, K_virtual, W_virtual] = predictYLZField( ...
    X_virtual, ...
    N_virtual, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

F_map = reshape(F_pred, size(xx));
Dispersion_map = reshape(F_dispersion, size(xx));

%% 7. 绘制 YLZ 核重构力分布图
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

% 绘制测点方向向量，便于观察各向异性方向
quiver( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    2.0 .* N_meas(:,1), ...
    2.0 .* N_meas(:,2), ...
    0, ...
    'k', ...
    'LineWidth', 1.2);

cb = colorbar;
ylabel(cb, 'Predicted force / N');

xlabel('x / mm');
ylabel('y / mm');

axis equal;
axis tight;
title('Force field reconstructed by YLZ anisotropic repulsive kernel');

%% 8. 绘制局部加权离散度
% 该图表示同一虚拟位置附近测点力值的加权差异。
% 数值越大，说明参与该位置预测的测点力值差异越明显。

figure;

contourf(xx, yy, Dispersion_map, 30, 'LineColor', 'none');
hold on;

scatter( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    50, ...
    'k', ...
    'filled');

cb = colorbar;
ylabel(cb, 'Local weighted dispersion / N');

xlabel('x / mm');
ylabel('y / mm');

axis equal;
axis tight;
title('Local weighted dispersion of YLZ kernel interpolation');

%% 9. 单个测点作为力源时的传播预测
% 单力源预测不能使用 K*F/K，否则核会被约掉，结果恒等于 F。
%
% 对单个测点，应使用：
%
%   F_single(x) = F_source * K_normalized(x)
%
% 这里选取中心测点作为示例，并将该测点对应的核归一化到最大值为 1。

source_index = 5;

K_single = K_virtual(:, source_index);
K_single_max = max(K_single);

if K_single_max > params.weight_tol
    K_single_normalized = K_single ./ K_single_max;
else
    K_single_normalized = zeros(size(K_single));
end

F_single = F_meas(source_index) .* K_single_normalized;
F_single_map = reshape(F_single, size(xx));

figure;

contourf(xx, yy, F_single_map, 30, 'LineColor', 'none');
hold on;

scatter( ...
    X_meas(source_index,1), ...
    X_meas(source_index,2), ...
    100, ...
    F_meas(source_index), ...
    'filled', ...
    'MarkerEdgeColor', 'k');

quiver( ...
    X_meas(source_index,1), ...
    X_meas(source_index,2), ...
    2.0 .* N_meas(source_index,1), ...
    2.0 .* N_meas(source_index,2), ...
    0, ...
    'k', ...
    'LineWidth', 1.5);

cb = colorbar;
ylabel(cb, 'Single-source predicted force / N');

xlabel('x / mm');
ylabel('y / mm');

axis equal;
axis tight;
title(sprintf( ...
    'Single-source YLZ force propagation: sensor %d', ...
    source_index));

%% 10. 输出所有整数坐标点的力数据
x_int = 0:1:20;
y_int = 0:1:20;

[xx_int, yy_int] = meshgrid(x_int, y_int);
X_int = [xx_int(:), yy_int(:)];

theta_int_deg = zeros(size(X_int, 1), 1);

N_int = [
    cosd(theta_int_deg), ...
    sind(theta_int_deg)
];

[F_int, F_int_dispersion] = predictYLZField( ...
    X_int, ...
    N_int, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

ResultTable = table( ...
    X_int(:,1), ...
    X_int(:,2), ...
    F_int, ...
    F_int_dispersion, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'PredictedForce_N', ...
        'LocalDispersion_N'} ...
);

disp('整数坐标点力数据：');
disp(ResultTable);

%% 11. 保存整数坐标点力数据到 Excel
output_file = 'integer_coordinate_force_result_YLZ.xlsx';

writetable(ResultTable, output_file);

fprintf('整数坐标点力数据已保存为 %s\n', output_file);

%% 12. 输出 YLZ 核参数
disp('YLZ anisotropic repulsive-kernel parameters:');
disp(params);

disp('说明：');
disp('1. 这些参数是人工设定的核参数，不是由当前程序自动学习得到。');
disp('2. LocalDispersion_N 是加权局部离散度，不是 GP 概率标准差。');
disp('3. 若所有方向向量相同，beta 的线性方向项为零。');
disp('4. 若需启用更强的方向不对称性，应为不同测点或虚拟点设置不同方向。');

%% ========================================================================
%  局部函数：使用 YLZ 各向异性排斥核进行力场预测
% ========================================================================

function [F_pred, F_dispersion, K, W] = predictYLZField( ...
    X_query, ...
    N_query, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params)

    % 检查输入尺寸
    if size(X_query,2) ~= 2 || size(X_meas,2) ~= 2
        error('X_query 和 X_meas 必须是两列二维坐标。');
    end

    if size(N_query,2) ~= 2 || size(N_meas,2) ~= 2
        error('N_query 和 N_meas 必须是两列二维方向向量。');
    end

    if size(X_query,1) ~= size(N_query,1)
        error('X_query 与 N_query 的行数必须一致。');
    end

    if size(X_meas,1) ~= size(N_meas,1)
        error('X_meas 与 N_meas 的行数必须一致。');
    end

    if size(X_meas,1) ~= numel(F_meas)
        error('测点数量必须与 F_meas 元素数量一致。');
    end

    % 将方向向量归一化，避免方向向量长度影响方向函数
    N_query = normalizeRows(N_query);
    N_meas = normalizeRows(N_meas);

    num_query = size(X_query, 1);
    num_meas = size(X_meas, 1);

    % K(q,j)：第 j 个真实测点对第 q 个虚拟点的 YLZ 核值
    K = zeros(num_query, num_meas);

    % D(q,j)：第 j 个真实测点到第 q 个虚拟点的距离
    D = zeros(num_query, num_meas);

    for j = 1:num_meas

        %% A. 几何距离与单位位移方向
        h = X_query - X_meas(j, :);

        r = sqrt(sum(h.^2, 2));
        D(:,j) = r;

        % 使用安全距离避免 r = 0 时除零
        r_safe = max(r, params.r_eps);

        % r_hat 指向：真实测点 j -> 当前虚拟点
        r_hat = bsxfun(@rdivide, h, r_safe);

        %% B. YLZ 径向余弦衰减部分
        u_radial = zeros(num_query, 1);

        near_region = r <= params.r_min;
        transition_region = ...
            r > params.r_min & r < params.r_c;

        % 在 r <= r_min 内保持最大势值，避免中心处奇异
        u_radial(near_region) = 2.0 .* params.epsilon0;

        % 在 r_min < r < r_c 内按 cos^(2*zeta) 衰减
        radial_argument = ...
            pi .* (r(transition_region) - params.r_min) ./ ...
            (2.0 .* (params.r_c - params.r_min));

        u_radial(transition_region) = ...
            2.0 .* params.epsilon0 .* ...
            cos(radial_argument).^(2.0 .* params.zeta);

        % 当 r >= r_c 时，u_radial 保持为 0

        %% C. YLZ 方向函数 a(r_hat, n_i, n_j)
        n_i = N_meas(j, :);
        n_j = N_query;

        % 二维叉乘只保留 z 分量：
        % [a_x,a_y,0] x [b_x,b_y,0] = [0,0,a_x*b_y-a_y*b_x]
        cross_i = ...
            n_i(1) .* r_hat(:,2) - ...
            n_i(2) .* r_hat(:,1);

        cross_j = ...
            n_j(:,1) .* r_hat(:,2) - ...
            n_j(:,2) .* r_hat(:,1);

        % (n_i - n_j)·r_hat
        direction_dot = ...
            (n_i(1) - n_j(:,1)) .* r_hat(:,1) + ...
            (n_i(2) - n_j(:,2)) .* r_hat(:,2);

        % a = (n_i x r_hat)*(n_j x r_hat)
        %     + beta*(n_i-n_j)·r_hat - beta^2
        a_value = ...
            cross_i .* cross_j + ...
            params.beta .* direction_dot - ...
            params.beta.^2;

        %% D. YLZ 方向调制函数 phi
        phi = 1.0 + params.mu .* (a_value - 1.0);

        % 在 r = 0 时，传播方向没有定义。
        % 这里令源点处 phi = 1，使核在源点保持最大且有限。
        coincident = r < params.coincide_tol;
        phi(coincident) = 1.0;

        % 作为插值权重时要求核非负。
        % 若需要严格保留原始有符号势函数，可删除此行，
        % 但可能出现负权重和接近零的权重和。
        phi = max(phi, 0.0);

        %% E. 完整 YLZ 各向异性排斥核
        K(:,j) = u_radial .* phi;
    end

    %% F. 对核权重进行归一化
    weight_sum = sum(K, 2);

    W = zeros(size(K));

    valid = weight_sum > params.weight_tol;

    W(valid,:) = bsxfun( ...
        @rdivide, ...
        K(valid,:), ...
        weight_sum(valid));

    %% G. 截止半径外无有效测点时，使用最近邻作为备用
    invalid_indices = find(~valid);

    if ~isempty(invalid_indices)
        [~, nearest_index] = min(D(invalid_indices,:), [], 2);

        linear_index = sub2ind( ...
            size(W), ...
            invalid_indices, ...
            nearest_index);

        W(linear_index) = 1.0;
    end

    %% H. 预测力值
    F_pred = W * F_meas;

    %% I. 若预测点与真实测点重合，强制返回真实测量值
    % 这样可保证测点位置处的结果严格等于对应测量值。
    [minimum_distance, nearest_measurement] = min(D, [], 2);

    exact_match = minimum_distance < params.coincide_tol;
    exact_indices = find(exact_match);

    if ~isempty(exact_indices)
        W(exact_indices,:) = 0.0;

        exact_linear_index = sub2ind( ...
            size(W), ...
            exact_indices, ...
            nearest_measurement(exact_indices));

        W(exact_linear_index) = 1.0;

        F_pred(exact_indices) = ...
            F_meas(nearest_measurement(exact_indices));
    end

    %% J. 计算加权局部离散度
    % 注意：该值不是概率意义上的标准差。
    force_difference = bsxfun( ...
        @minus, ...
        F_meas.', ...
        F_pred);

    F_dispersion = sqrt( ...
        sum(W .* force_difference.^2, 2));

    F_dispersion(exact_match) = 0.0;
end

%% ========================================================================
%  局部函数：逐行归一化二维方向向量
% ========================================================================

function N_unit = normalizeRows(N)

    row_norm = sqrt(sum(N.^2, 2));

    if any(row_norm < eps)
        error('方向向量不能为零向量。');
    end

    N_unit = bsxfun(@rdivide, N, row_norm);
end
