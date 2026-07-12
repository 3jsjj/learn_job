clc;
clear;
close all;

%% 1. 真实霍尔测点三维坐标 (mm)
% 每一行均为 [x, y, z]。
X_meas = [
     0,  0, 0;
    10,  0, 0;
    20,  0, 0;
     0, 10, 0;
    10, 10, 0;
    20, 10, 0;
     0, 20, 0;
    10, 20, 0;
    20, 20, 0
];

%% 2. 霍尔芯片计算出的真实力 (N)
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

%% 3. YLZ 三维各向异性排斥核参数
% 径向部分：
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
%   a = (n_t x r_hat)·(n_v x r_hat)
%       + beta*(n_t-n_v)·r_hat
%       - beta^2
%
% 其中：
%   n_t   = 真实测点方向向量
%   n_v   = 当前虚拟点方向向量
%   r_hat = 从真实测点指向当前虚拟点的单位向量
%
% 完整 YLZ 核：
%
%   K_qj = u_R(r_qj)*phi_qj
%
% 当前虚拟点的预测力：
%
%               sum_{j in neighbors(q)} K_qj*F_j
%   F_pred(q) = ------------------------------------
%                    sum_{j in neighbors(q)} K_qj
%
% neighbors(q) 只包含三维直线距离小于 r_c 的真实测点。

params.epsilon0 = 1.0;     % 径向势公共幅值
params.r_min    = 2.0;     % 近场平台半径 (mm)
params.r_c      = 12.0;    % 三维直线距离截止半径 rcut (mm)
params.zeta     = 2.0;     % 径向衰减陡峭程度
params.mu       = 0.8;     % 各向异性强度
params.beta     = 0.2;     % 方向偏置参数

params.r_eps        = 1e-12;  % 防止 r = 0 时计算 r_hat 出现除零
params.weight_tol   = 1e-12;  % 判断核权重和是否有效
params.coincide_tol = 1e-10;  % 判断虚拟点是否与真实测点重合

% 如果某个虚拟点在 r_c 内没有真实邻居：
% true  -> 使用三维最近邻真实测点作为备用；
% false -> 返回 NaN，严格遵守 r_c 截断。
params.use_nearest_fallback = true;

if params.r_c <= params.r_min
    error('参数错误：r_c 必须大于 r_min。');
end

%% 4. 根据三维坐标推导真实测点方向
% 单独一个坐标点不能唯一确定方向，因此需要指定方向参考点 O。
%
% 当前采用径向方向模型：
%
%   n_t = (X_t - O)/||X_t - O||
%
% 这里把方向参考点放在测点几何中心下方。
% direction_origin_depth 应根据实际结构中的等效方向中心进行调整。

direction_origin_depth = 5.0;  % mm

direction_origin = [
    mean(X_meas(:,1)), ...
    mean(X_meas(:,2)), ...
    min(X_meas(:,3)) - direction_origin_depth
];

[N_meas, azimuth_meas_deg, elevation_meas_deg] = ...
    directionsFromCoordinates( ...
        X_meas, ...
        direction_origin, ...
        params.r_eps);

%% 5. 创建三维连续预测网格
x = linspace(0, 20, 41);
y = linspace(0, 20, 41);
z = linspace(0, 10, 21);

[xx, yy, zz] = meshgrid(x, y, z);

% 每行是一个虚拟点 [x, y, z]
X_virtual = [ xx(:), yy(:), zz(:)];

%% 6. 根据虚拟点坐标推导虚拟点方向
% 使用与真实测点相同的方向参考点：
%
%   n_v = (X_v - O)/||X_v - O||

[N_virtual, azimuth_virtual_deg, elevation_virtual_deg] = ...
    directionsFromCoordinates( ...
        X_virtual, ...
        direction_origin, ...
        params.r_eps);

%% 7. 使用“虚拟点 -> neighbors -> 真实点”的双层循环预测
% 外层循环：
%   逐个处理虚拟点 q。
%
% 内层循环：
%   只处理三维直线距离 r < r_c 的真实邻居 j。
%
% 函数内部会明确计算：
%   n_v、n_t、h_tv、r、r_hat、u_R、a、phi、K_qj，
% 然后再对当前虚拟点的全部邻居进行归一化加权。

[
    F_pred, ...
    F_dispersion, ...
    K_virtual, ...
    W_virtual, ...
    NeighborList, ...
    GeometricNeighborCount, ...
    ActiveKernelCount ...
] = predictYLZField3DNeighbors( ...
    X_virtual, ...
    N_virtual, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

F_volume = reshape(F_pred, size(xx));
Dispersion_volume = reshape(F_dispersion, size(xx));

GeometricNeighborCount_volume = ...
    reshape(GeometricNeighborCount, size(xx));

ActiveKernelCount_volume = ...
    reshape(ActiveKernelCount, size(xx));

%% 8. 绘制三维预测力场切片
x_slices = mean(x);
y_slices = mean(y);
z_slices = unique([ ...
    z(1), ...
    z(round((numel(z)+1)/2)), ...
    z(end) ...
]);

figure;

slice( ...
    xx, yy, zz, F_volume, ...
    x_slices, y_slices, z_slices);

shading interp;
hold on;

scatter3( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    90, ...
    F_meas, ...
    'filled', ...
    'MarkerEdgeColor', 'k');

direction_scale = 2.0;

quiver3( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    direction_scale .* N_meas(:,1), ...
    direction_scale .* N_meas(:,2), ...
    direction_scale .* N_meas(:,3), ...
    0, ...
    'k', ...
    'LineWidth', 1.2);

cb = colorbar;
ylabel(cb, 'Predicted force / N');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title('3D YLZ force field: virtual-point neighbor-loop prediction');

%% 9. 绘制局部加权离散度
figure;

slice( ...
    xx, yy, zz, Dispersion_volume, ...
    x_slices, y_slices, z_slices);

shading interp;
hold on;

scatter3( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    55, ...
    'k', ...
    'filled');

cb = colorbar;
ylabel(cb, 'Local weighted dispersion / N');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title('Local weighted dispersion');

%% 10. 绘制 r_c 内的几何邻居数量
% 该数量只由三维直线距离 r < r_c 决定，
% 不考虑 phi 截断后核值是否变为零。

figure;

slice( ...
    xx, yy, zz, GeometricNeighborCount_volume, ...
    x_slices, y_slices, z_slices);

shading flat;
hold on;

scatter3( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    55, ...
    'k', ...
    'filled');

cb = colorbar;
ylabel(cb, 'Geometric neighbor count');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title('Number of real sensors inside r_c');

%% 11. 绘制最终非零核数量
% 即满足 r < r_c 且 K_qj > weight_tol 的真实测点数量。

figure;

slice( ...
    xx, yy, zz, ActiveKernelCount_volume, ...
    x_slices, y_slices, z_slices);

shading flat;
hold on;

scatter3( ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    55, ...
    'k', ...
    'filled');

cb = colorbar;
ylabel(cb, 'Active kernel count');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title('Number of nonzero YLZ kernels');

%% 12. 单个真实测点作为三维力源时的传播
% 单力源传播不能使用 K*F/K，否则核会被约掉。
%
% 应使用：
%
%   F_single = F_source*K_normalized

source_index = 5;

K_single = K_virtual(:, source_index);
K_single_max = max(K_single);

if K_single_max > params.weight_tol
    K_single_normalized = K_single ./ K_single_max;
else
    K_single_normalized = zeros(size(K_single));
end

F_single = F_meas(source_index) .* K_single_normalized;
F_single_volume = reshape(F_single, size(xx));

figure;

slice( ...
    xx, yy, zz, F_single_volume, ...
    x_slices, y_slices, z_slices);

shading interp;
hold on;

scatter3( ...
    X_meas(source_index,1), ...
    X_meas(source_index,2), ...
    X_meas(source_index,3), ...
    110, ...
    F_meas(source_index), ...
    'filled', ...
    'MarkerEdgeColor', 'k');

quiver3( ...
    X_meas(source_index,1), ...
    X_meas(source_index,2), ...
    X_meas(source_index,3), ...
    direction_scale .* N_meas(source_index,1), ...
    direction_scale .* N_meas(source_index,2), ...
    direction_scale .* N_meas(source_index,3), ...
    0, ...
    'k', ...
    'LineWidth', 1.5);

cb = colorbar;
ylabel(cb, 'Single-source predicted force / N');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title(sprintf( ...
    'Single-source YLZ propagation: sensor %d', ...
    source_index));

%% 13. 输出所有整数三维坐标点的预测结果
x_int = 0:1:20;
y_int = 0:1:20;
z_int = 0:1:10;

[xx_int, yy_int, zz_int] = meshgrid( ...
    x_int, ...
    y_int, ...
    z_int);

X_int = [
    xx_int(:), ...
    yy_int(:), ...
    zz_int(:)
];

[N_int, azimuth_int_deg, elevation_int_deg] = ...
    directionsFromCoordinates( ...
        X_int, ...
        direction_origin, ...
        params.r_eps);

[
    F_int, ...
    F_int_dispersion, ...
    K_int, ...
    W_int, ...
    NeighborList_int, ...
    GeometricNeighborCount_int, ...
    ActiveKernelCount_int ...
] = predictYLZField3DNeighbors( ...
    X_int, ...
    N_int, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

ResultTable = table( ...
    X_int(:,1), ...
    X_int(:,2), ...
    X_int(:,3), ...
    F_int, ...
    F_int_dispersion, ...
    GeometricNeighborCount_int, ...
    ActiveKernelCount_int, ...
    azimuth_int_deg, ...
    elevation_int_deg, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'z_mm', ...
        'PredictedForce_N', ...
        'LocalDispersion_N', ...
        'NeighborCountInsideRcut', ...
        'ActiveKernelCount', ...
        'Azimuth_deg', ...
        'Elevation_deg'} ...
);

disp('整数三维坐标点力数据：');
disp(ResultTable);

%% 14. 保存整数三维坐标点力数据到 Excel
output_file = ...
    'integer_3D_coordinate_force_result_YLZ_neighbor_loop.xlsx';

writetable(ResultTable, output_file);

fprintf( ...
    '整数三维坐标点力数据已保存为 %s\n', ...
    output_file);

%% 15. 输出真实测点方向与核参数
DirectionTable = table( ...
    (1:size(X_meas,1)).', ...
    X_meas(:,1), ...
    X_meas(:,2), ...
    X_meas(:,3), ...
    N_meas(:,1), ...
    N_meas(:,2), ...
    N_meas(:,3), ...
    azimuth_meas_deg, ...
    elevation_meas_deg, ...
    'VariableNames', { ...
        'SensorIndex', ...
        'x_mm', ...
        'y_mm', ...
        'z_mm', ...
        'n_tx', ...
        'n_ty', ...
        'n_tz', ...
        'Azimuth_deg', ...
        'Elevation_deg'} ...
);

disp('真实测点方向：');
disp(DirectionTable);

disp('YLZ kernel parameters:');
disp(params);

fprintf( ...
    '方向参考点 direction_origin = [%.6f, %.6f, %.6f] mm\n', ...
    direction_origin(1), ...
    direction_origin(2), ...
    direction_origin(3));

disp('说明：');
disp('1. 外层循环逐个处理虚拟点 q。');
disp('2. 每个虚拟点先用三维直线距离筛选 r < r_c 的真实邻居。');
disp('3. 内层循环逐个计算 n_v、n_t、r_hat、a、phi 和 K_qj。');
disp('4. 最终只对当前 neighbors 的核预测结果进行归一化加权。');
disp('5. LocalDispersion_N 不是 GP 概率标准差。');

%% ========================================================================
%  局部函数：
%  对每一个虚拟点先找 r_c 内的真实 neighbors，
%  再逐个计算 YLZ 核，最后进行归一化加权预测
% ========================================================================

function [ F_pred, F_dispersion, K, W, NeighborList, GeometricNeighborCount, ...
    ActiveKernelCount ] = predictYLZField3DNeighbors( X_query, N_query, ...
    X_meas, N_meas, F_meas, params)

    %% A. 输入检查
    if size(X_query,2) ~= 3 || size(X_meas,2) ~= 3
        error('X_query 和 X_meas 必须是三列坐标 [x,y,z]。');
    end

    if size(N_query,2) ~= 3 || size(N_meas,2) ~= 3
        error('N_query 和 N_meas 必须是三列方向向量 [nx,ny,nz]。');
    end

    if size(X_query,1) ~= size(N_query,1)
        error('每个虚拟点必须有一个对应的方向向量。');
    end

    if size(X_meas,1) ~= size(N_meas,1)
        error('每个真实测点必须有一个对应的方向向量。');
    end

    if size(X_meas,1) ~= numel(F_meas)
        error('真实测点数量必须与真实力数量一致。');
    end

    if params.r_c <= params.r_min
        error('参数错误：r_c 必须大于 r_min。');
    end

    %% B. 预处理
    N_query = normalizeRows(N_query);
    N_meas = normalizeRows(N_meas);
    F_meas = F_meas(:);

    num_query = size(X_query, 1);
    num_meas = size(X_meas, 1);

    F_pred = nan(num_query, 1);
    F_dispersion = nan(num_query, 1);

    % 完整矩阵仍然保留，便于后续绘图和检查。
    % 非 neighbor 的元素始终为 0。
    K = zeros(num_query, num_meas);
    W = zeros(num_query, num_meas);

    % 每个虚拟点对应的真实邻居编号
    NeighborList = cell(num_query, 1);

    % r < r_c 的几何邻居数
    GeometricNeighborCount = zeros(num_query, 1);

    % 最终 K_qj > weight_tol 的有效核数量
    ActiveKernelCount = zeros(num_query, 1);

    %% C. 外层循环：逐个虚拟点进行预测
    for q = 1:num_query

        % ---------------------------------------------------------------
        % C1. 取出当前虚拟点的位置和方向
        % ---------------------------------------------------------------
        x_v = X_query(q, :);
        n_v = N_query(q, :);

        % ---------------------------------------------------------------
        % C2. 计算当前虚拟点到所有真实测点的三维直线距离
        % ---------------------------------------------------------------
        % delta_all(j,:) = x_v - x_t(j,:)
        %
        % 即从真实测点 j 指向当前虚拟点 q 的三维位移。
        delta_all = bsxfun( ...
            @minus, ...
            x_v, ...
            X_meas);

        distance_all = sqrt( ...
            sum(delta_all.^2, 2));

        % ---------------------------------------------------------------
        % C3. 根据三维直线距离筛选 r_c 内的真实 neighbors
        % ---------------------------------------------------------------
        neighbors = find( ...
            distance_all < params.r_c);

        NeighborList{q} = neighbors(:).';

        GeometricNeighborCount(q) = ...
            numel(neighbors);

        % ---------------------------------------------------------------
        % C4. 当前虚拟点没有任何 r_c 内邻居时的处理
        % ---------------------------------------------------------------
        if isempty(neighbors)

            if params.use_nearest_fallback

                % 使用三维直线距离最近的真实测点作为备用。
                [~, nearest_j] = min(distance_all);

                F_pred(q) = F_meas(nearest_j);
                F_dispersion(q) = 0.0;

                W(q, nearest_j) = 1.0;
                ActiveKernelCount(q) = 0;

            else

                % 严格执行 r_c 截断：没有 neighbor 就不预测。
                F_pred(q) = NaN;
                F_dispersion(q) = NaN;

            end

            continue;
        end

        % ---------------------------------------------------------------
        % C5. 若当前虚拟点与某个真实测点重合，直接返回真实力
        % ---------------------------------------------------------------
        exact_neighbors = neighbors( ...
            distance_all(neighbors) < ...
            params.coincide_tol);

        if ~isempty(exact_neighbors)

            % 如果存在重复坐标，只取第一个重合真实测点。
            exact_j = exact_neighbors(1);

            K(q, exact_j) = ...
                2.0 .* params.epsilon0;

            W(q, exact_j) = 1.0;

            F_pred(q) = F_meas(exact_j);
            F_dispersion(q) = 0.0;

            ActiveKernelCount(q) = 1;

            continue;
        end

        % 当前虚拟点所有邻居的核值。
        local_kernel = zeros(numel(neighbors), 1);

        % 分子：
        %   sum_j K_qj*F_j
        weighted_force_sum = 0.0;

        % 分母：
        %   sum_j K_qj
        kernel_sum = 0.0;

        % ---------------------------------------------------------------
        % C6. 内层循环：逐个处理当前虚拟点的真实 neighbor
        % ---------------------------------------------------------------
        for local_index = 1:numel(neighbors)

            % 当前真实邻居的全局编号
            j = neighbors(local_index);

            % -----------------------------------------------------------
            % C6.1 取出真实点位置、真实力和真实点方向 n_t
            % -----------------------------------------------------------
            x_t = X_meas(j, :);
            F_t = F_meas(j);
            n_t = N_meas(j, :);

            % -----------------------------------------------------------
            % C6.2 从真实点指向虚拟点的位移、距离和 r_hat
            % -----------------------------------------------------------
            % h_tv 的方向是：
            %   true point t -> virtual point v
            h_tv = x_v - x_t;

            % 三维直线距离
            r = sqrt( ...
                h_tv(1)^2 + ...
                h_tv(2)^2 + ...
                h_tv(3)^2);

            % 防止除以 0
            r_safe = max(r, params.r_eps);

            % 从真实点指向虚拟点的三维单位向量
            r_hat = h_tv ./ r_safe;

            % -----------------------------------------------------------
            % C6.3 计算 YLZ 径向排斥部分 u_R(r)
            % -----------------------------------------------------------
            if r <= params.r_min

                % 近场平台
                u_radial = ...
                    2.0 .* params.epsilon0;

            elseif r < params.r_c

                % r_min < r < r_c 内平滑余弦衰减
                radial_argument = ...
                    pi .* (r - params.r_min) ./ ...
                    (2.0 .* ...
                    (params.r_c - params.r_min));

                u_radial = ...
                    2.0 .* params.epsilon0 .* ...
                    cos(radial_argument).^( ...
                    2.0 .* params.zeta);

            else

                % 理论上不会进入该分支，因为 neighbors 已经过滤。
                u_radial = 0.0;

            end

            % -----------------------------------------------------------
            % C6.4 计算 YLZ 方向函数中的叉乘项
            % -----------------------------------------------------------
            % n_t x r_hat
            cross_t = cross(n_t, r_hat);

            % n_v x r_hat
            cross_v = cross(n_v, r_hat);

            % (n_t x r_hat)·(n_v x r_hat)
            cross_dot = dot( ...
                cross_t, ...
                cross_v);

            % -----------------------------------------------------------
            % C6.5 计算方向差沿 r_hat 的投影
            % -----------------------------------------------------------
            % (n_t - n_v)·r_hat
            direction_dot = dot( ...
                n_t - n_v, ...
                r_hat);

            % -----------------------------------------------------------
            % C6.6 计算 YLZ 中间方向量 a
            % -----------------------------------------------------------
            a_value = ...
                cross_dot + ...
                params.beta .* direction_dot - ...
                params.beta.^2;

            % -----------------------------------------------------------
            % C6.7 计算方向调制函数 phi
            % -----------------------------------------------------------
            phi = ...
                1.0 + ...
                params.mu .* (a_value - 1.0);

            % 作为归一化插值权重时要求非负。
            phi = max(phi, 0.0);

            % -----------------------------------------------------------
            % C6.8 计算完整 YLZ 核 K_qj
            % -----------------------------------------------------------
            K_qj = ...
                u_radial .* phi;

            % 保存当前虚拟点 q 与真实点 j 的核值
            K(q, j) = K_qj;
            local_kernel(local_index) = K_qj;

            % -----------------------------------------------------------
            % C6.9 真实力 F_t 通过核 K_qj 向虚拟点预测
            % -----------------------------------------------------------
            % 当前真实测点对虚拟点产生的加权力贡献：
            %
            %   contribution_qj = K_qj*F_t
            contribution_qj = ...
                K_qj .* F_t;

            % 累积分子和分母
            weighted_force_sum = ...
                weighted_force_sum + ...
                contribution_qj;

            kernel_sum = ...
                kernel_sum + ...
                K_qj;
        end

        % ---------------------------------------------------------------
        % C7. 对当前虚拟点所有 neighbor 的预测进行归一化
        % ---------------------------------------------------------------
        ActiveKernelCount(q) = sum( ...
            local_kernel > params.weight_tol);

        if kernel_sum > params.weight_tol

            % 归一化加权预测
            F_pred(q) = ...
                weighted_force_sum ./ kernel_sum;

            % 保存归一化权重
            W(q, neighbors) = ...
                (local_kernel ./ kernel_sum).';

            % -----------------------------------------------------------
            % C8. 计算当前虚拟点的局部加权离散度
            % -----------------------------------------------------------
            local_weights = ...
                W(q, neighbors).';

            local_force_difference = ...
                F_meas(neighbors) - ...
                F_pred(q);

            F_dispersion(q) = sqrt( ...
                sum( ...
                    local_weights .* ...
                    local_force_difference.^2));
        else

            % 虽然存在几何 neighbor，但其 phi 全部被截断为 0。
            if params.use_nearest_fallback

                [~, local_nearest_index] = min( ...
                    distance_all(neighbors));

                nearest_j = neighbors( ...
                    local_nearest_index);

                F_pred(q) = F_meas(nearest_j);
                F_dispersion(q) = 0.0;

                W(q, nearest_j) = 1.0;

            else

                F_pred(q) = NaN;
                F_dispersion(q) = NaN;

            end
        end
    end
end

%% ========================================================================
%  局部函数：
%  根据三维坐标和公共方向参考点计算单位方向、方位角和仰角
% ========================================================================

function [
    N_unit, ...
    azimuth_deg, ...
    elevation_deg ...
] = directionsFromCoordinates( ...
    X, ...
    origin, ...
    zero_tol)

    if size(X,2) ~= 3
        error('X 必须是三列三维坐标 [x,y,z]。');
    end

    if ~isvector(origin) || numel(origin) ~= 3
        error('origin 必须是三维参考点 [x0,y0,z0]。');
    end

    origin = reshape(origin, 1, 3);

    % 从公共参考点指向每个坐标点的向量
    V = bsxfun( ...
        @minus, ...
        X, ...
        origin);

    radius = sqrt( ...
        sum(V.^2, 2));

    if any(radius < zero_tol)

        bad_index = find( ...
            radius < zero_tol, ...
            1, ...
            'first');

        error( ...
            ['方向参考点与第 %d 个坐标点重合，' ...
             '无法确定方向，请移动 direction_origin。'], ...
            bad_index);
    end

    % 单位方向向量
    N_unit = bsxfun( ...
        @rdivide, ...
        V, ...
        radius);

    % 方位角：x-y 平面投影相对 +x 轴的角度
    azimuth_deg = atan2d( ...
        V(:,2), ...
        V(:,1));

    % 仰角：向量相对 x-y 平面的角度
    horizontal_radius = hypot( ...
        V(:,1), ...
        V(:,2));

    elevation_deg = atan2d( ...
        V(:,3), ...
        horizontal_radius);
end

%% ========================================================================
%  局部函数：逐行归一化三维方向向量
% ========================================================================

function N_unit = normalizeRows(N)

    row_norm = sqrt( ...
        sum(N.^2, 2));

    if any(row_norm < eps)
        error('方向向量不能为零向量。');
    end

    N_unit = bsxfun( ...
        @rdivide, ...
        N, ...
        row_norm);
end
