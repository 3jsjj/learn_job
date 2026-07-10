clc;
clear;
close all;

%% 1. 真实霍尔测点三维坐标 (mm)
% 每一行均为 [x, y, z]。
%
% 当前示例中的 9 个测点都位于 z = 0 平面。
% 程序可以在三维空间中预测，但由于训练测点只分布在一个平面上，
% z 方向的变化完全由人工指定的 YLZ 核决定，而不是由多层实测数据学习得到。
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
% 三维方向部分：
%
%   phi = 1 + mu*(a - 1)
%
%   a = (n_source x r_hat)·(n_query x r_hat)
%       + beta*(n_source-n_query)·r_hat
%       - beta^2
%
% 完整核：
%
%   K = u_R(r)*phi
%
% 为了将势函数作为非负插值权重，程序执行：
%
%   phi = max(phi, 0)

params.epsilon0 = 1.0;     % 径向势整体幅值；归一化插值时公共幅值会约掉
params.r_min    = 2.0;     % 近场平台半径 (mm)
params.r_c      = 12.0;    % 截止半径 (mm)，必须大于 r_min
params.zeta     = 2.0;     % 径向衰减陡峭程度
params.mu       = 0.8;     % 各向异性强度
params.beta     = 0.2;     % 方向偏置参数

params.r_eps        = 1e-12;  % 避免 r = 0 时单位方向向量除零
params.weight_tol   = 1e-12;  % 判断核权重和是否有效
params.coincide_tol = 1e-10;  % 判断查询点是否与真实测点重合

if params.r_c <= params.r_min
    error('参数错误：r_c 必须大于 r_min。');
end

%% 4. 设置真实测点的三维方向向量
% 三维单位方向向量使用方位角 azimuth 和仰角 elevation 表示：
%
%   n_x = cos(elevation)*cos(azimuth)
%   n_y = cos(elevation)*sin(azimuth)
%   n_z = sin(elevation)
%
% 角度单位均为 degree。
%
% azimuth = 0°，elevation = 0°  对应 +x 方向；
% azimuth = 90°，elevation = 0° 对应 +y 方向；
% elevation = 90°               对应 +z 方向。
%
% 当前示例假设所有真实测点的方向均沿 +z。
% 如果每个测点具有不同方向，应逐行设置下面两个角度向量。

num_meas = size(X_meas, 1);

azimuth_meas_deg   = zeros(num_meas, 1);
elevation_meas_deg = 90 .* ones(num_meas, 1);

N_meas = directionFromAzEl( ...
    azimuth_meas_deg, ...
    elevation_meas_deg);

%% 5. 创建三维连续预测网格
% x、y、z 的范围和采样数量可根据实际结构修改。
%
% 当前测点均位于 z = 0。这里示例预测 z = 0~10 mm 的空间。
% 如果需要预测测点上下两侧，可改成：
%
%   z = linspace(-10, 10, 41);

x = linspace(0, 20, 41);
y = linspace(0, 20, 41);
z = linspace(0, 10, 21);

% meshgrid 生成尺寸为：
%   length(y) × length(x) × length(z)
[xx, yy, zz] = meshgrid(x, y, z);

% 每行是一个三维虚拟查询点 [x, y, z]
X_virtual = [xx(:), yy(:), zz(:)];

%% 6. 设置三维虚拟点的方向向量
% 当前假设所有虚拟点的局部方向也。
%
% 若虚拟点方向随空间位置变化，需要分别构造：
%   azimuth_virtual_deg
%   elevation_virtual_deg

num_virtual = size(X_virtual, 1);
% 这里的方位角都市为 0°，仰角为 90°，即所有虚拟点方向均沿 +z。
azimuth_virtual_deg   = zeros(num_virtual, 1);
elevation_virtual_deg = 90 .* ones(num_virtual, 1);

N_virtual = directionFromAzEl( ...
    azimuth_virtual_deg, ...
    elevation_virtual_deg);

%% 7. 使用三维 YLZ 各向异性排斥核预测连续力场
% 多测点重构采用归一化核加权：
%
%                sum_j K(q,j)*F_meas(j)
%   F_pred(q) = -------------------------
%                      sum_j K(q,j)
%
% F_dispersion 是局部加权离散度，不是概率意义上的标准差。

[F_pred, F_dispersion, K_virtual, W_virtual] = predictYLZField3D( ...
    X_virtual, ...
    N_virtual, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

% 将列向量恢复为三维体数据
F_volume = reshape(F_pred, size(xx));
Dispersion_volume = reshape(F_dispersion, size(xx));

% 每个虚拟点处具有非零核值的真实测点数量
ActiveCount = sum(K_virtual > params.weight_tol, 2);
ActiveCount_volume = reshape(ActiveCount, size(xx));

%% 8. 绘制三维预测力场切片
% 三维标量场不能直接使用二维 contourf 完整显示。
% 这里使用 slice 绘制 x、y、z 三组切片。

x_slices = mean(x);
y_slices = mean(y);
z_slices = unique([z(1), z(round((numel(z)+1)/2)), z(end)]);

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

title('3D force field reconstructed by YLZ anisotropic repulsive kernel');

%% 9. 绘制三维局部加权离散度切片
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

title('3D local weighted dispersion of YLZ interpolation');

%% 10. 绘制有效测点数量切片
% 如果某处 ActiveCount = 0，说明该位置位于所有测点的截止半径之外，
% 程序会使用最近邻测点作为备用预测。

figure;

slice( ...
    xx, yy, zz, ActiveCount_volume, ...
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
ylabel(cb, 'Number of active sensors');

xlabel('x / mm');
ylabel('y / mm');
zlabel('z / mm');

axis equal;
axis tight;
grid on;
view(40, 28);

title('Number of active YLZ kernels in 3D space');

%% 11. 单个测点作为三维力源时的传播预测
% 单力源传播不能使用 K*F/K，否则 K 会被约掉，结果恒等于 F。
%
% 应使用：
%
%   F_single(x,y,z) = F_source*K_normalized(x,y,z)

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
    '3D single-source YLZ propagation: sensor %d', ...
    source_index));

%% 12. 输出所有整数三维坐标点的力数据
% 这里的整数 z 范围与连续预测空间保持一致。
% 根据实际需求修改 z_int。

x_int = 0:1:20;
y_int = 0:1:20;
z_int = 0:1:10;

[xx_int, yy_int, zz_int] = meshgrid(x_int, y_int, z_int);

X_int = [
    xx_int(:), ...
    yy_int(:), ...
    zz_int(:)
];

num_int = size(X_int, 1);

azimuth_int_deg   = zeros(num_int, 1);
elevation_int_deg = 90 .* ones(num_int, 1);

N_int = directionFromAzEl( ...
    azimuth_int_deg, ...
    elevation_int_deg);

[F_int, F_int_dispersion, K_int] = predictYLZField3D( ...
    X_int, ...
    N_int, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params);

ActiveSensorCount_int = sum(K_int > params.weight_tol, 2);

ResultTable = table( ...
    X_int(:,1), ...
    X_int(:,2), ...
    X_int(:,3), ...
    F_int, ...
    F_int_dispersion, ...
    ActiveSensorCount_int, ...
    'VariableNames', { ...
        'x_mm', ...
        'y_mm', ...
        'z_mm', ...
        'PredictedForce_N', ...
        'LocalDispersion_N', ...
        'ActiveSensorCount'} ...
);

disp('整数三维坐标点力数据：');
disp(ResultTable);

%% 13. 保存整数三维坐标点力数据到 Excel
output_file = 'integer_3D_coordinate_force_result_YLZ.xlsx';

writetable(ResultTable, output_file);

fprintf('整数三维坐标点力数据已保存为 %s\n', output_file);

%% 14. 输出 YLZ 核参数与方向设置
disp('YLZ 3D anisotropic repulsive-kernel parameters:');
disp(params);

disp('真实测点方向向量 N_meas：');
disp(N_meas);

disp('说明：');
disp('1. 所有坐标、距离和方向计算均已改为三维。');
disp('2. LocalDispersion_N 是加权局部离散度，不是 GP 概率标准差。');
disp('3. 若所有源点和查询点方向相同，beta 的线性方向项为零。');
disp('4. 当前测点均位于 z=0，因此 z 方向预测是核模型外推。');
disp('5. 若要从数据标定 z 方向衰减，应增加不同 z 高度的真实测点。');

%% ========================================================================
%  局部函数：使用三维 YLZ 各向异性排斥核进行力场预测
% ========================================================================

function [F_pred, F_dispersion, K, W] = predictYLZField3D( ...
    X_query, ...
    N_query, ...
    X_meas, ...
    N_meas, ...
    F_meas, ...
    params)

    %% A. 检查三维输入尺寸
    if size(X_query,2) ~= 3 || size(X_meas,2) ~= 3
        error('X_query 和 X_meas 必须是三列三维坐标 [x,y,z]。');
    end

    if size(N_query,2) ~= 3 || size(N_meas,2) ~= 3
        error('N_query 和 N_meas 必须是三列三维方向向量 [nx,ny,nz]。');
    end

    if size(X_query,1) ~= size(N_query,1)
        error('X_query 与 N_query 的行数必须一致。');
    end

    if size(X_meas,1) ~= size(N_meas,1)
        error('X_meas 与 N_meas 的行数必须一致。');
    end

    if size(X_meas,1) ~= numel(F_meas)
        error('真实测点数量必须与 F_meas 元素数量一致。');
    end

    if params.r_c <= params.r_min
        error('参数错误：r_c 必须大于 r_min。');
    end

    %% B. 将三维方向向量逐行归一化
    N_query = normalizeRows(N_query);
    N_meas = normalizeRows(N_meas);

    F_meas = F_meas(:);

    num_query = size(X_query, 1);
    num_meas = size(X_meas, 1);

    % K(q,j)：真实测点 j 对查询点 q 的三维 YLZ 核值
    K = zeros(num_query, num_meas);

    % D(q,j)：真实测点 j 到查询点 q 的三维欧氏距离
    D = zeros(num_query, num_meas);

    %% C. 逐个真实测点计算其对全部查询点的三维核值
    for j = 1:num_meas

        % ---------------------------------------------------------------
        % C1. 三维位移、距离和单位传播方向
        % ---------------------------------------------------------------
        % h(q,:) = [x_q-x_j, y_q-y_j, z_q-z_j]
        h = bsxfun(@minus, X_query, X_meas(j,:));

        % 三维欧氏距离：
        % r = sqrt(dx^2 + dy^2 + dz^2)
        r = sqrt(sum(h.^2, 2));
        D(:,j) = r;

        % 防止 r = 0 时单位方向向量除零
        r_safe = max(r, params.r_eps);

        % 从真实测点 j 指向查询点 q 的单位向量
        r_hat = bsxfun(@rdivide, h, r_safe);

        % ---------------------------------------------------------------
        % C2. YLZ 径向有限支撑余弦核
        % ---------------------------------------------------------------
        u_radial = zeros(num_query, 1);

        near_region = r <= params.r_min;

        transition_region = ...
            r > params.r_min & ...
            r < params.r_c;

        % 近场平台
        u_radial(near_region) = ...
            2.0 .* params.epsilon0;

        % r_min 到 r_c 之间平滑衰减
        radial_argument = ...
            pi .* (r(transition_region) - params.r_min) ./ ...
            (2.0 .* (params.r_c - params.r_min));

        u_radial(transition_region) = ...
            2.0 .* params.epsilon0 .* ...
            cos(radial_argument).^(2.0 .* params.zeta);

        % r >= r_c 的位置保持为 0

        % ---------------------------------------------------------------
        % C3. 三维 YLZ 方向函数
        % ---------------------------------------------------------------
        n_source = N_meas(j,:);
        n_source_all = repmat(n_source, num_query, 1);
        n_query = N_query;

        % 三维叉乘：
        % cross_source(q,:) = n_source x r_hat(q,:)
        % cross_query(q,:)  = n_query(q,:) x r_hat(q,:)
        cross_source = cross(n_source_all, r_hat, 2);
        cross_query  = cross(n_query,      r_hat, 2);

        % 两个三维叉乘向量的点积
        cross_dot = sum( ...
            cross_source .* cross_query, ...
            2);

        % 方向差沿传播方向的投影：
        % (n_source - n_query)·r_hat
        direction_difference = ...
            n_source_all - n_query;

        direction_dot = sum( ...
            direction_difference .* r_hat, ...
            2);

        % YLZ 方向中间量
        a_value = ...
            cross_dot + ...
            params.beta .* direction_dot - ...
            params.beta.^2;

        % ---------------------------------------------------------------
        % C4. 方向调制函数 phi
        % ---------------------------------------------------------------
        phi = ...
            1.0 + ...
            params.mu .* (a_value - 1.0);

        % 当查询点与源点重合时，r_hat 无确定方向。
        % 令 phi = 1，使源点处核值保持有限且不受方向项压低。
        coincident = r < params.coincide_tol;
        phi(coincident) = 1.0;

        % 作为插值权重时，要求核值非负。
        phi = max(phi, 0.0);

        % ---------------------------------------------------------------
        % C5. 完整三维 YLZ 各向异性排斥核
        % ---------------------------------------------------------------
        K(:,j) = u_radial .* phi;
    end

    %% D. 对每个查询点的核权重归一化
    weight_sum = sum(K, 2);

    W = zeros(size(K));

    valid = weight_sum > params.weight_tol;

    W(valid,:) = bsxfun( ...
        @rdivide, ...
        K(valid,:), ...
        weight_sum(valid));

    %% E. 所有核均为零时使用三维最近邻作为备用
    invalid_indices = find(~valid);

    if ~isempty(invalid_indices)

        [~, nearest_index] = min( ...
            D(invalid_indices,:), ...
            [], ...
            2);

        linear_index = sub2ind( ...
            size(W), ...
            invalid_indices, ...
            nearest_index);

        W(linear_index) = 1.0;
    end

    %% F. 计算预测力
    F_pred = W * F_meas;

    %% G. 查询点与真实测点重合时强制返回真实测量值
    [minimum_distance, nearest_measurement] = min(D, [], 2);

    exact_match = ...
        minimum_distance < params.coincide_tol;

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

    %% H. 计算局部加权离散度
    % D_F(q) = sqrt(sum_j W(q,j)*(F_j-F_pred(q))^2)
    %
    % 该值反映参与同一预测的真实测点力值差异，
    % 不是概率意义上的预测标准差。

    force_difference = bsxfun( ...
        @minus, ...
        F_meas.', ...
        F_pred);

    F_dispersion = sqrt( ...
        sum( ...
            W .* force_difference.^2, ...
            2));

    F_dispersion(exact_match) = 0.0;
end

%% ========================================================================
%  局部函数：将方位角和仰角转换为三维单位方向向量
% ========================================================================

function N = directionFromAzEl(azimuth_deg, elevation_deg)

    azimuth_deg = azimuth_deg(:);
    elevation_deg = elevation_deg(:);

    if numel(azimuth_deg) ~= numel(elevation_deg)
        error('方位角和仰角的元素数量必须一致。');
    end

    N = [
        cosd(elevation_deg) .* cosd(azimuth_deg), ...
        cosd(elevation_deg) .* sind(azimuth_deg), ...
        sind(elevation_deg)
    ];

    N = normalizeRows(N);
end

%% ========================================================================
%  局部函数：逐行归一化方向向量
% ========================================================================

function N_unit = normalizeRows(N)

    row_norm = sqrt(sum(N.^2, 2));

    if any(row_norm < eps)
        error('方向向量不能为零向量。');
    end

    N_unit = bsxfun(@rdivide, N, row_norm);
end
