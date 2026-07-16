% =========================================================================
% 曲面压力逆有限元重构（简洁版）
%
% 目标：
%   1) 输入真实测点位移和真实测点力/压力；
%   2) 在 virtual_coords 指定的曲面节点上反演标量压力 p；
%   3) 根据重构压力计算虚拟节点位移；
%
% 曲面压力不是独立的 Fx、Fy、Fz。
% 每个虚拟节点只有一个未知压力 p_i，并转换为等效节点力：
%
%       f_i = pressure_sign * p_i * A_i * n_i
%
% 其中：
%   n_i：曲面节点单位法向；
%   A_i：节点分摊面积；
%   pressure_sign = -1 表示正压力沿外法向反方向作用。
%
% 输入文件类型与原程序一致：
%   .inp  : Abaqus 节点和单元
%   .mtx  : 5 列刚度矩阵 [NodeI, DOFI, NodeJ, DOFJ, Value]
%   .csv  : 测点数据、固定节点、可选虚拟点坐标
%
% 测点 CSV 支持两种格式：
%   格式 A：NodeID, U1, U2, U3, F1, F2, F3
%   格式 B：NodeID, U1, U2, U3, Pressure
%
% 单位：
%   若刚度为 N/mm、坐标为 mm，则输出压力为 N/mm^2，即 MPa。
%
% 本版本适用于当前上传的壳单元 .mtx：
%   每节点 6 自由度，压力仅装配到 U1/U2/U3 平移自由度。
% =========================================================================

clear;
clc;

%% 1. 用户设置

inp_filepath      = 'element_nodes_2.inp';
mtx_filepath      = 'get_matrix-3_STIF2.mtx';
measured_filepath = 'Abaqus_Nodal_U_and_Pressure.csv';
% 固定节点坐标文件：
%   推荐格式 A：[OriginalNodeID, X, Y, Z]
%   也支持格式 B：[X, Y, Z]
% 坐标必须与 inp 中 node_coords 使用同一坐标系和同一长度单位。
fixed_filepath = 'fixed_nodes_with_coords.csv';

% 固定坐标与 inp 节点坐标的匹配容差。
% 若坐标单位为 mm，通常可先尝试 1e-6～1e-4。
fixed_coordinate_tolerance = 1e-6;

% 可选文件：每行 [X, Y, Z]，必须对应目标曲面的有限元节点。
% 如果文件不存在，程序默认使用模型全部外表面节点。
virtual_coords_filepath = 'virtual_coords.csv';

% 当前 .mtx 已确认是壳单元矩阵：每个节点 6 个自由度
% 1-3 为平移 U1/U2/U3，4-6 为转动 UR1/UR2/UR3
num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];

% 若 fixed_nodes.csv 表示“完全固支边缘”，约束全部 6 个自由度。
% 如果边界不是完全固支，请改成实际受约束的自由度编号。
fixed_dof_components = 1:6;

% +1：正压力沿计算出的表面法向
% -1：正压力沿表面法向的反方向，通常用于外部压力压向结构
pressure_sign = -1;

% 虚拟坐标到有限元节点的绝对匹配容差
coordinate_tolerance = 1e-6;

% 相对 Tikhonov 正则化参数
% 可尝试 1e-8、1e-6、1e-4、1e-2
lambda_relative = 1e-6;

% 若测点已有已知力/压力，通常不再把测点压力作为未知量
exclude_measured_nodes_from_virtual = true;


%% 2. 读取节点、单元及模型外表面三角形

[node_ids, node_coords, exterior_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath);

fprintf('读取节点数量：%d\n', numel(node_ids));
fprintf('识别外表面三角形数量：%d\n', size(exterior_triangles, 1));


%% 3. 读取 virtual_coords 并映射到有限元节点

if isfile(virtual_coords_filepath)
    virtual_coords = readmatrix(virtual_coords_filepath);
    virtual_coords = virtual_coords( ...
        all(isfinite(virtual_coords), 2), :);

    if size(virtual_coords, 2) == 2
        virtual_coords(:, 3) = 0;
    end

    if size(virtual_coords, 2) ~= 3
        error('virtual_coords.csv 必须包含 2 列或 3 列坐标。');
    end
else
    % 没有单独的虚拟坐标文件时，默认使用全部外表面节点
    virtual_node_rows_default = unique(exterior_triangles(:));
    virtual_coords = node_coords(virtual_node_rows_default, :);

    warning(['未找到 virtual_coords.csv，', ...
             '当前使用模型全部外表面节点作为虚拟压力节点。']);
end

[virtual_node_rows, virtual_match_distance] = ...
    match_coordinates_to_nodes(virtual_coords, node_coords);

valid_coordinate_match = ...
    virtual_match_distance <= coordinate_tolerance;

if any(~valid_coordinate_match)
    warning('%d 个虚拟坐标未匹配到有限元节点，已剔除。', ...
        nnz(~valid_coordinate_match));
end

virtual_node_rows = virtual_node_rows(valid_coordinate_match);
virtual_coords    = virtual_coords(valid_coordinate_match, :);

% 去除多个坐标映射到同一个有限元节点的情况
[virtual_node_rows, unique_virtual_index] = ...
    unique(virtual_node_rows, 'stable');
virtual_coords = virtual_coords(unique_virtual_index, :);

if isempty(virtual_node_rows)
    error('没有有效的虚拟曲面节点。');
end


%% 4. 读取真实测点位移和力/压力

measured_data = readmatrix(measured_filepath);

if size(measured_data, 2) < 5
    error(['测点文件至少应为：', ...
           'NodeID, U1, U2, U3, Pressure，', ...
           '或者 NodeID, U1, U2, U3, F1, F2, F3。']);
end

id_measured_raw = measured_data(:, 1);
U_measured_raw  = measured_data(:, 2:4);

if size(measured_data, 2) >= 7
    measured_input_type = 'force';
    measured_load_raw = measured_data(:, 5:7);

    valid_measured_data = ...
        isfinite(id_measured_raw) & ...
        all(isfinite(U_measured_raw), 2) & ...
        all(isfinite(measured_load_raw), 2);
else
    measured_input_type = 'pressure';
    measured_load_raw = measured_data(:, 5);

    valid_measured_data = ...
        isfinite(id_measured_raw) & ...
        all(isfinite(U_measured_raw), 2) & ...
        isfinite(measured_load_raw);
end

id_measured_raw = id_measured_raw(valid_measured_data);
U_measured_raw  = U_measured_raw(valid_measured_data, :);
measured_load_raw = measured_load_raw(valid_measured_data, :);

[measured_id_exists, measured_node_rows_raw] = ...
    ismember(id_measured_raw, node_ids);

if any(~measured_id_exists)
    warning('%d 个测点 NodeID 不在 inp 节点中，已剔除。', ...
        nnz(~measured_id_exists));
end

id_measured_raw    = id_measured_raw(measured_id_exists);
U_measured_raw     = U_measured_raw(measured_id_exists, :);
measured_load_raw  = measured_load_raw(measured_id_exists, :);
measured_node_rows = measured_node_rows_raw(measured_id_exists);

if isempty(measured_node_rows)
    error('没有有效的真实测点。');
end


%% 5. 从目标曲面网格计算节点法向和分摊面积

% 目标曲面应由 virtual_coords 所覆盖。
% 测点若位于同一曲面，也纳入曲面节点集合。
target_surface_node_rows = unique( ...
    [virtual_node_rows; measured_node_rows]);

is_target_triangle = all( ...
    ismember(exterior_triangles, target_surface_node_rows), 2);

target_surface_triangles = ...
    exterior_triangles(is_target_triangle, :);

if isempty(target_surface_triangles)
    error([ ...
        '无法由 virtual_coords 形成曲面单元。', newline, ...
        'virtual_coords 应包含目标曲面的有限元表面节点，', ...
        '而不是任意空间散点。']);
end

[node_normals, node_areas] = compute_nodal_surface_geometry( ...
    node_coords, target_surface_triangles);

% 检查虚拟节点是否真正属于该曲面
valid_virtual_surface = ...
    node_areas(virtual_node_rows) > 0 & ...
    vecnorm(node_normals(virtual_node_rows, :), 2, 2) > 0;

if any(~valid_virtual_surface)
    warning('%d 个虚拟节点没有有效曲面面积或法向，已剔除。', ...
        nnz(~valid_virtual_surface));
end

virtual_node_rows = virtual_node_rows(valid_virtual_surface);
virtual_coords    = virtual_coords(valid_virtual_surface, :);


%% 6. 读取并组装 Abaqus 刚度矩阵

% .mtx 在 MATLAB 中不是 readmatrix 自动识别的标准扩展名，
% 因此明确指定它是文本文件。若当前 MATLAB 版本不支持该写法，
% 则回退到 load（适用于纯数字五列表格）。
try
    mtx_data = readmatrix(mtx_filepath, 'FileType', 'text');
catch
    mtx_data = load(mtx_filepath);
end

mtx_data = mtx_data(all(isfinite(mtx_data), 2), :);

if size(mtx_data, 2) < 5
    error('刚度矩阵文件必须为 5 列格式。');
end

node_i = mtx_data(:, 1);
dof_i  = mtx_data(:, 2);
node_j = mtx_data(:, 3);
dof_j  = mtx_data(:, 4);
values = mtx_data(:, 5);

matrix_dof_components = unique([dof_i; dof_j])';

fprintf('mtx 中出现的节点自由度编号：');
fprintf('%d ', matrix_dof_components);
fprintf('\n');

if any(~ismember(matrix_dof_components, 1:num_dofs_per_node))
    error(['mtx 中存在超出 num_dofs_per_node 的自由度编号。', ...
           '请检查每节点自由度数设置。']);
end

if ~all(ismember(translation_dofs, matrix_dof_components))
    error('mtx 中缺少 U1/U2/U3 平移自由度。');
end

row_idx = num_dofs_per_node .* (node_i - 1) + dof_i;
col_idx = num_dofs_per_node .* (node_j - 1) + dof_j;

max_dof = max([row_idx; col_idx]);

K_input = sparse( ...
    row_idx, col_idx, values, max_dof, max_dof);

% 判断输入是半矩阵还是完整矩阵
has_upper = nnz(triu(K_input, 1)) > 0;
has_lower = nnz(tril(K_input, -1)) > 0;

if xor(has_upper, has_lower)
    K_raw = K_input + K_input' ...
        - spdiags(diag(K_input), 0, max_dof, max_dof);
else
    symmetry_error = norm(K_input - K_input', 'fro') / ...
        max(norm(K_input, 'fro'), eps);

    if symmetry_error > 1e-8
        warning('输入刚度矩阵不是严格对称矩阵，当前进行对称化。');
    end

    K_raw = (K_input + K_input') / 2;
end


%% 7. 根据固定节点坐标映射本地 NodeID，并施加边界条件

fixed_data = readmatrix(fixed_filepath);
fixed_data = fixed_data(all(isfinite(fixed_data), 2), :);

if size(fixed_data, 2) >= 4
    % [OriginalNodeID, X, Y, Z]
    fixed_original_ids = fixed_data(:, 1);
    fixed_coords = fixed_data(:, 2:4);

elseif size(fixed_data, 2) == 3
    % [X, Y, Z]
    fixed_original_ids = nan(size(fixed_data, 1), 1);
    fixed_coords = fixed_data(:, 1:3);

elseif size(fixed_data, 2) == 2 && size(node_coords, 2) == 3
    % 二维坐标输入时补 Z=0
    fixed_original_ids = nan(size(fixed_data, 1), 1);
    fixed_coords = [fixed_data(:, 1:2), zeros(size(fixed_data, 1), 1)];

else
    error([ ...
        '固定节点文件必须包含坐标：', newline, ...
        '[X,Y,Z] 或 [OriginalNodeID,X,Y,Z]。']);
end

% 根据坐标找到 inp 中最近的有限元节点行号
[fixed_node_rows, fixed_match_distance] = ...
    match_coordinates_to_nodes(fixed_coords, node_coords);

fixed_match_ok = ...
    fixed_match_distance <= fixed_coordinate_tolerance;

fprintf('\n');
fprintf('============== 固定节点坐标映射 ==============\n');
fprintf('输入固定坐标数量：%d\n', size(fixed_coords, 1));
fprintf('成功匹配数量：%d\n', nnz(fixed_match_ok));
fprintf('匹配失败数量：%d\n', nnz(~fixed_match_ok));
fprintf('最大最近节点距离：%.12e\n', max(fixed_match_distance));
fprintf('平均最近节点距离：%.12e\n', mean(fixed_match_distance));
fprintf('允许匹配容差：%.12e\n', fixed_coordinate_tolerance);
fprintf('==============================================\n');

if any(~fixed_match_ok)
    failed_rows = find(~fixed_match_ok);
    print_count = min(10, numel(failed_rows));

    fprintf('前 %d 个匹配失败固定点：\n', print_count);

    for k = 1:print_count
        r = failed_rows(k);

        fprintf([ ...
            '  行 %d，输入坐标=[%.9g, %.9g, %.9g]，', ...
            '最近距离=%.6e\n'], ...
            r, ...
            fixed_coords(r, 1), ...
            fixed_coords(r, 2), ...
            fixed_coords(r, 3), ...
            fixed_match_distance(r));
    end

    error([ ...
        '部分固定节点坐标无法匹配到 inp 节点。', newline, ...
        '请检查坐标系、单位、实例平移/旋转以及匹配容差。']);
end

% 获取 inp 使用的本地节点 ID
fixed_nodes_local = node_ids(fixed_node_rows);

% 防止多个输入坐标匹配到同一个节点
[fixed_nodes_local, unique_fixed_index] = ...
    unique(fixed_nodes_local, 'stable');

fixed_node_rows = fixed_node_rows(unique_fixed_index);
fixed_coords = fixed_coords(unique_fixed_index, :);
fixed_original_ids = fixed_original_ids(unique_fixed_index);

% 确认这些本地 NodeID 确实存在于 mtx 中
mtx_node_ids = unique([node_i; node_j]);
fixed_in_mtx = ismember(fixed_nodes_local, mtx_node_ids);

fprintf('去重后的固定节点数量：%d\n', ...
    numel(fixed_nodes_local));
fprintf('本地 NodeID 能在 mtx 中找到：%d / %d\n', ...
    nnz(fixed_in_mtx), numel(fixed_nodes_local));

if any(~fixed_in_mtx)
    bad_local_ids = fixed_nodes_local(~fixed_in_mtx);

    error([ ...
        '%d 个坐标匹配得到的 inp NodeID 不在 mtx 中。', newline, ...
        '这说明 inp 与 mtx 可能不是同一网格或使用不同编号。'], ...
        numel(bad_local_ids));
end

% 壳单元完全固支：每个固定节点删除 1～6 自由度
constrained_dofs_matrix = node_dof_numbers( ...
    fixed_nodes_local, ...
    num_dofs_per_node, ...
    fixed_dof_components);

constrained_dofs = constrained_dofs_matrix(:);

[nonzero_row, nonzero_col] = find(K_raw);
existing_dofs = unique([nonzero_row; nonzero_col]);

constrained_dofs_existing = intersect( ...
    constrained_dofs, existing_dofs);

active_dofs = setdiff( ...
    existing_dofs, constrained_dofs_existing);

K_global = K_raw(active_dofs, active_dofs);

fprintf('\n');
fprintf('刚度矩阵原始非零自由度：%d\n', ...
    numel(existing_dofs));
fprintf('固定节点数量：%d\n', ...
    numel(fixed_nodes_local));
fprintf('理论应删除的约束自由度：%d\n', ...
    numel(fixed_nodes_local) * numel(fixed_dof_components));
fprintf('实际删除的约束自由度：%d\n', ...
    numel(constrained_dofs_existing));
fprintf('活动自由度数量：%d\n', ...
    size(K_global, 1));

% 绘制固定节点映射结果，进行肉眼验证
figure( ...
    'Name', '固定节点坐标映射检查', ...
    'Color', 'w', ...
    'Position', [100, 100, 900, 700]);

trisurf( ...
    exterior_triangles, ...
    node_coords(:, 1), ...
    node_coords(:, 2), ...
    node_coords(:, 3), ...
    'FaceColor', [0.75, 0.80, 0.90], ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', [0.65, 0.65, 0.65]);

hold on;

scatter3( ...
    node_coords(fixed_node_rows, 1), ...
    node_coords(fixed_node_rows, 2), ...
    node_coords(fixed_node_rows, 3), ...
    35, ...
    'r', ...
    'filled');

axis equal;
grid on;
view(45, 30);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('红色节点应与 Abaqus 中的固定边缘完全一致');
legend('模型表面', '坐标映射后的固定节点');
hold off;


%% 8. 过滤不在活动自由度中的测点

measured_global_dofs = node_dof_numbers( ...
    id_measured_raw, ...
    num_dofs_per_node, ...
    translation_dofs);

valid_measured_active = all( ...
    ismember(measured_global_dofs, active_dofs), 2);

if any(~valid_measured_active)
    warning('%d 个测点属于固定自由度或无刚度自由度，已剔除。', ...
        nnz(~valid_measured_active));
end

id_measured       = id_measured_raw(valid_measured_active);
measured_node_rows = measured_node_rows(valid_measured_active);
U_measured_mat    = U_measured_raw(valid_measured_active, :);
measured_load     = measured_load_raw(valid_measured_active, :);

measured_global_dofs = measured_global_dofs( ...
    valid_measured_active, :);

[~, dof_measured_matrix] = ismember( ...
    measured_global_dofs, active_dofs);

dof_measured = reshape(dof_measured_matrix', [], 1);
U_measured   = reshape(U_measured_mat', [], 1);

if isempty(dof_measured)
    error('过滤后没有有效测点自由度。');
end


%% 9. 将测点已知力或已知压力转换为节点力

if strcmp(measured_input_type, 'force')
    % CSV 已直接给出 F1、F2、F3
    F_measured_mat = measured_load;
else
    % CSV 给出标量压力，转换为沿局部曲面法向的等效节点力
    measured_normals = node_normals(measured_node_rows, :);
    measured_areas   = node_areas(measured_node_rows);

    if any(measured_areas <= 0)
        error('部分压力测点没有有效的曲面分摊面积。');
    end

    measured_force_per_pressure = bsxfun( ...
        @times, measured_normals, pressure_sign .* measured_areas);

    F_measured_mat = bsxfun( ...
        @times, measured_force_per_pressure, measured_load);
end

F_measured = reshape(F_measured_mat', [], 1);

N_total_dof = size(K_global, 1);

F_known = sparse( ...
    dof_measured, ...
    ones(size(dof_measured)), ...
    F_measured, ...
    N_total_dof, ...
    1);


%% 10. 过滤虚拟压力节点并构造压力到节点力的映射 G

id_virtual = node_ids(virtual_node_rows);

virtual_global_dofs = node_dof_numbers( ...
    id_virtual, ...
    num_dofs_per_node, ...
    translation_dofs);

% 将各种过滤原因分别统计，避免把“与测点重复”误认为“固定”
virtual_has_stiffness = all( ...
    ismember(virtual_global_dofs, existing_dofs), 2);

virtual_is_fixed = ismember( ...
    id_virtual, fixed_nodes_local);

virtual_is_active = all( ...
    ismember(virtual_global_dofs, active_dofs), 2);

virtual_is_measured = ismember( ...
    id_virtual, id_measured);

fprintf('\n');
fprintf('============== 虚拟节点诊断 ==============\n');
fprintf('初始虚拟表面节点：%d\n', numel(id_virtual));
fprintf('列在 fixed_nodes.csv 中：%d\n', ...
    nnz(virtual_is_fixed));
fprintf('U1/U2/U3 未完整出现在 mtx 中：%d\n', ...
    nnz(~virtual_has_stiffness));
fprintf('与真实测点 NodeID 重复：%d\n', ...
    nnz(virtual_is_measured));
fprintf('具有活动自由度：%d\n', ...
    nnz(virtual_is_active));

if exclude_measured_nodes_from_virtual
    valid_virtual_active = ...
        virtual_is_active & ~virtual_is_measured;
else
    valid_virtual_active = virtual_is_active;
end

fprintf('最终保留的虚拟节点：%d\n', ...
    nnz(valid_virtual_active));
fprintf('==========================================\n');

if any(virtual_is_fixed)
    warning('%d 个虚拟节点属于固定边界，已剔除。', ...
        nnz(virtual_is_fixed));
end

if any(~virtual_has_stiffness)
    warning('%d 个虚拟节点没有完整的刚度自由度，已剔除。', ...
        nnz(~virtual_has_stiffness));
end

if exclude_measured_nodes_from_virtual && any(virtual_is_measured)
    warning('%d 个虚拟节点与真实测点重复，已剔除。', ...
        nnz(virtual_is_measured));
end

id_virtual        = id_virtual(valid_virtual_active);
virtual_node_rows = virtual_node_rows(valid_virtual_active);
virtual_coords    = virtual_coords(valid_virtual_active, :);
virtual_global_dofs = virtual_global_dofs(valid_virtual_active, :);

[~, dof_virtual_matrix] = ismember( ...
    virtual_global_dofs, active_dofs);

N_virtual_nodes = numel(id_virtual);

if N_virtual_nodes == 0
    if exclude_measured_nodes_from_virtual && ...
            all(virtual_is_measured | ~virtual_is_active)

        error([ ...
            '过滤后没有可反演的虚拟压力节点。', newline, ...
            '最常见原因是测点 CSV 包含了全部表面节点，', ...
            '同时 exclude_measured_nodes_from_virtual=true。', ...
            newline, ...
            '请查看上方“虚拟节点诊断”的分类计数。']);
    else
        error([ ...
            '过滤后没有可反演的虚拟压力节点。', newline, ...
            '请查看上方“固定节点诊断”和“虚拟节点诊断”。']);
    end
end

virtual_normals = node_normals(virtual_node_rows, :);
virtual_areas   = node_areas(virtual_node_rows);

% 每单位压力在节点三个方向产生的等效力
% 尺寸：[N_virtual_nodes, 3]
force_per_unit_pressure = bsxfun( ...
    @times, virtual_normals, pressure_sign .* virtual_areas);

% 构造稀疏映射：
% F_pressure_global = G * pressure_virtual
G_rows = reshape(dof_virtual_matrix', [], 1);
G_cols = repelem((1:N_virtual_nodes)', 3);
G_vals = reshape(force_per_unit_pressure', [], 1);

G = sparse( ...
    G_rows, G_cols, G_vals, ...
    N_total_dof, N_virtual_nodes);


%% 11. 构造“虚拟压力 -> 测点位移”灵敏度矩阵

% E_measured 从全局活动位移中提取测点自由度
N_measured_dof = numel(dof_measured);

E_measured = sparse( ...
    dof_measured, ...
    (1:N_measured_dof)', ...
    ones(N_measured_dof, 1), ...
    N_total_dof, ...
    N_measured_dof);

% 伴随形式：
% H = E_measured' * inv(K_global) * G
% 这样右端项数量由测点自由度数决定
Z = K_global' \ E_measured;
H = Z' * G;


%% 12. 扣除测点已知载荷产生的基准位移

U_from_known_force = K_global \ F_known;
U_measured_baseline = U_from_known_force(dof_measured);

Delta_U = U_measured - U_measured_baseline;


%% 13. Tikhonov 正则化反演曲面压力

% 求解：
% min ||H*p - Delta_U||^2 + lambda*||p||^2
%
% 使用对偶形式，避免构造巨大的 H'*H

H_scale = norm(H, 'fro')^2 / max(N_virtual_nodes, 1);
lambda_absolute = lambda_relative * max(H_scale, eps);

dual_matrix = ...
    H * H' + lambda_absolute * speye(N_measured_dof);

pressure_virtual = ...
    H' * (dual_matrix \ Delta_U);


%% 14. 根据重构压力计算完整位移和虚拟点位移

F_pressure_global = G * pressure_virtual;
F_total_global    = F_known + F_pressure_global;

U_reconstructed_global = K_global \ F_total_global;

dof_virtual = reshape(dof_virtual_matrix', [], 1);

U_virtual_vector = ...
    U_reconstructed_global(dof_virtual);

U_virtual_3D = reshape( ...
    U_virtual_vector, 3, [])';

% 压力产生的三维等效节点力，仅用于检查
F_virtual_equivalent = bsxfun( ...
    @times, force_per_unit_pressure, pressure_virtual);


%% 15. 误差诊断

U_measured_predicted = ...
    U_reconstructed_global(dof_measured);

relative_displacement_error = norm( ...
    U_measured_predicted - U_measured) / ...
    max(norm(U_measured), eps);

equilibrium_residual = norm( ...
    K_global * U_reconstructed_global - F_total_global) / ...
    max(norm(F_total_global), eps);

fprintf('\n');
fprintf('================ 重构诊断 ================\n');
fprintf('有效测点数量：%d\n', numel(id_measured));
fprintf('虚拟压力节点数量：%d\n', N_virtual_nodes);
fprintf('测点位移相对拟合误差：%.6e\n', ...
    relative_displacement_error);
fprintf('全局平衡相对残差：%.6e\n', ...
    equilibrium_residual);
fprintf('最小重构压力：%.6e\n', min(pressure_virtual));
fprintf('最大重构压力：%.6e\n', max(pressure_virtual));
fprintf('==========================================\n');


%% 16. 输出 CSV

result_table = table( ...
    id_virtual(:), ...
    virtual_coords(:, 1), ...
    virtual_coords(:, 2), ...
    virtual_coords(:, 3), ...
    virtual_normals(:, 1), ...
    virtual_normals(:, 2), ...
    virtual_normals(:, 3), ...
    virtual_areas(:), ...
    pressure_virtual(:), ...
    U_virtual_3D(:, 1), ...
    U_virtual_3D(:, 2), ...
    U_virtual_3D(:, 3), ...
    F_virtual_equivalent(:, 1), ...
    F_virtual_equivalent(:, 2), ...
    F_virtual_equivalent(:, 3), ...
    'VariableNames', { ...
        'NodeID', ...
        'X', 'Y', 'Z', ...
        'NormalX', 'NormalY', 'NormalZ', ...
        'NodalArea', ...
        'Pressure', ...
        'U1_Reconstructed', ...
        'U2_Reconstructed', ...
        'U3_Reconstructed', ...
        'EquivalentForceX', ...
        'EquivalentForceY', ...
        'EquivalentForceZ'});

writetable( ...
    result_table, ...
    'Curved_Surface_Pressure_Reconstruction.csv');

disp('结果已保存：Curved_Surface_Pressure_Reconstruction.csv');


%% 17. 曲面压力点云显示

figure( ...
    'Name', '曲面压力逆有限元重构', ...
    'Color', 'w', ...
    'Position', [100, 100, 900, 700]);

scatter3( ...
    virtual_coords(:, 1), ...
    virtual_coords(:, 2), ...
    virtual_coords(:, 3), ...
    45, pressure_virtual, 'filled');

axis equal;
grid on;
view(45, 30);

xlabel('X');
ylabel('Y');
zlabel('Z');
title('曲面虚拟节点压力重构');

cb = colorbar;
ylabel(cb, 'Pressure');

colormap jet;

hold on;

model_length = max( ...
    max(node_coords, [], 1) - min(node_coords, [], 1));

normal_display_length = 0.03 * model_length;

quiver3( ...
    virtual_coords(:, 1), ...
    virtual_coords(:, 2), ...
    virtual_coords(:, 3), ...
    virtual_normals(:, 1), ...
    virtual_normals(:, 2), ...
    virtual_normals(:, 3), ...
    normal_display_length, ...
    'Color', [0.2, 0.2, 0.2]);

hold off;


% =========================================================================
%                              局部函数
% =========================================================================

function dof_matrix = node_dof_numbers( ...
    node_id_list, dofs_per_node, components)
% 将节点 ID 和指定自由度分量转换为绝对自由度编号。
%
% 壳单元每节点 6 自由度时：
%   node_dof_numbers(5, 6, [1 2 3]) -> [25 26 27]
%   node_dof_numbers(5, 6, 1:6)     -> [25 26 27 28 29 30]

    if nargin < 3
        components = 1:dofs_per_node;
    end

    node_id_list = node_id_list(:);
    components = components(:)';

    if any(components < 1) || ...
            any(components > dofs_per_node) || ...
            any(components ~= round(components))
        error('自由度分量 components 设置无效。');
    end

    dof_matrix = bsxfun( ...
        @plus, ...
        dofs_per_node .* (node_id_list - 1), ...
        components);
end


function [matched_rows, distances] = ...
    match_coordinates_to_nodes(query_coords, node_coords)
% 将坐标匹配到最近的有限元节点。
% 先尝试完全匹配，再对未匹配坐标进行最近邻搜索。
% 不依赖 Statistics Toolbox。

    N_query = size(query_coords, 1);

    matched_rows = zeros(N_query, 1);
    distances    = inf(N_query, 1);

    [exact_match, exact_rows] = ...
        ismember(query_coords, node_coords, 'rows');

    matched_rows(exact_match) = exact_rows(exact_match);
    distances(exact_match) = 0;

    unmatched = find(~exact_match);

    for k = 1:numel(unmatched)
        q = unmatched(k);

        delta = bsxfun( ...
            @minus, node_coords, query_coords(q, :));

        squared_distance = sum(delta.^2, 2);

        [minimum_squared_distance, nearest_row] = ...
            min(squared_distance);

        matched_rows(q) = nearest_row;
        distances(q) = sqrt(minimum_squared_distance);
    end
end


function [node_normals, node_areas] = ...
    compute_nodal_surface_geometry(node_coords, triangles)
% 根据曲面三角形计算：
%   node_normals：面积加权节点单位法向
%   node_areas  ：每个三角形面积平均分给三个节点

    N_nodes = size(node_coords, 1);

    normal_accumulator = zeros(N_nodes, 3);
    node_areas = zeros(N_nodes, 1);

    for t = 1:size(triangles, 1)
        rows = triangles(t, :);

        x1 = node_coords(rows(1), :);
        x2 = node_coords(rows(2), :);
        x3 = node_coords(rows(3), :);

        area_vector = cross(x2 - x1, x3 - x1);
        twice_area  = norm(area_vector);

        if twice_area <= eps
            continue;
        end

        triangle_area = 0.5 * twice_area;

        % 每个节点获得三角形面积的 1/3
        node_areas(rows) = ...
            node_areas(rows) + triangle_area / 3;

        % area_vector / 6 =
        % 单位法向 * triangle_area / 3
        normal_contribution = area_vector / 6;

        normal_accumulator(rows, :) = ...
            normal_accumulator(rows, :) + ...
            repmat(normal_contribution, 3, 1);
    end

    node_normals = zeros(N_nodes, 3);

    normal_length = vecnorm( ...
        normal_accumulator, 2, 2);

    valid_normal = normal_length > 0;

    node_normals(valid_normal, :) = bsxfun( ...
        @rdivide, ...
        normal_accumulator(valid_normal, :), ...
        normal_length(valid_normal));
end


function [node_ids, node_coords, surface_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath)
% 读取常见 Abaqus 节点和单元，并提取表面三角形。
%
% 当前支持：
%   壳单元：S3、S3R、S4、S4R
%   实体单元：C3D4、C3D8、C3D8R
%
% 对实体单元：
%   自动删除两个实体共享的内部面，只保留外表面。
%
% surface_triangles 中保存的是 node_coords 的行号，而不是 NodeID。

    raw_text = fileread(inp_filepath);
    lines = splitlines(raw_text);

    node_ids = zeros(0, 1);
    node_coords = zeros(0, 3);

    element_types = cell(0, 1);
    element_connectivity = cell(0, 1);

    current_mode = '';
    current_element_type = '';
    expected_element_nodes = 0;
    element_buffer = [];

    for line_index = 1:numel(lines)
        line = strtrim(lines{line_index});

        if isempty(line) || startsWith(line, '**')
            continue;
        end

        if startsWith(line, '*')
            upper_line = upper(line);

            if startsWith(upper_line, '*NODE')
                current_mode = 'node';
                current_element_type = '';
                expected_element_nodes = 0;
                element_buffer = [];

            elseif startsWith(upper_line, '*ELEMENT')
                current_mode = 'element';

                token = regexp( ...
                    upper_line, ...
                    'TYPE\s*=\s*([^,\s]+)', ...
                    'tokens', 'once');

                if isempty(token)
                    current_element_type = '';
                    expected_element_nodes = 0;
                else
                    current_element_type = token{1};
                    expected_element_nodes = ...
                        nodes_per_supported_element( ...
                            current_element_type);
                end

                element_buffer = [];
            else
                current_mode = '';
                current_element_type = '';
                expected_element_nodes = 0;
                element_buffer = [];
            end

            continue;
        end

        values = sscanf( ...
            strrep(line, ',', ' '), '%f')';

        if isempty(values)
            continue;
        end

        if strcmp(current_mode, 'node')
            if numel(values) < 3
                continue;
            end

            node_ids(end + 1, 1) = values(1);

            xyz = zeros(1, 3);
            number_of_coordinates = min(numel(values) - 1, 3);
            xyz(1:number_of_coordinates) = ...
                values(2:1 + number_of_coordinates);

            node_coords(end + 1, :) = xyz;

        elseif strcmp(current_mode, 'element') && ...
                expected_element_nodes > 0

            element_buffer = [element_buffer, values]; %#ok<AGROW>

            values_per_element = ...
                expected_element_nodes + 1;

            while numel(element_buffer) >= values_per_element
                one_element = ...
                    element_buffer(1:values_per_element);

                element_buffer(1:values_per_element) = [];

                element_types{end + 1, 1} = ...
                    current_element_type;

                element_connectivity{end + 1, 1} = ...
                    one_element(2:end);
            end
        end
    end

    if isempty(node_ids)
        error('未能从 inp 文件读取节点。');
    end

    if isempty(element_connectivity)
        error(['未读取到受支持的单元。', newline, ...
               '当前支持 S3/S4/C3D4/C3D8。']);
    end

    id_to_row = containers.Map( ...
        'KeyType', 'double', ...
        'ValueType', 'double');

    for n = 1:numel(node_ids)
        id_to_row(node_ids(n)) = n;
    end

    shell_faces = cell(0, 1);

    face_count = containers.Map( ...
        'KeyType', 'char', ...
        'ValueType', 'double');

    face_oriented = containers.Map( ...
        'KeyType', 'char', ...
        'ValueType', 'any');

    for element_index = 1:numel(element_connectivity)
        element_type = element_types{element_index};
        connectivity = element_connectivity{element_index};

        if startsWith(element_type, 'S3')
            shell_faces{end + 1, 1} = connectivity(1:3);

        elseif startsWith(element_type, 'S4')
            shell_faces{end + 1, 1} = connectivity(1:4);

        elseif startsWith(element_type, 'C3D4')
            local_faces = { ...
                [1, 2, 3], ...
                [1, 4, 2], ...
                [2, 4, 3], ...
                [3, 4, 1]};

            add_solid_faces( ...
                connectivity, local_faces, ...
                node_coords, id_to_row, ...
                face_count, face_oriented);

        elseif startsWith(element_type, 'C3D8')
            local_faces = { ...
                [1, 2, 3, 4], ...
                [5, 8, 7, 6], ...
                [1, 5, 6, 2], ...
                [2, 6, 7, 3], ...
                [3, 7, 8, 4], ...
                [4, 8, 5, 1]};

            add_solid_faces( ...
                connectivity, local_faces, ...
                node_coords, id_to_row, ...
                face_count, face_oriented);
        end
    end

    surface_faces = shell_faces;

    solid_face_keys = keys(face_count);

    for key_index = 1:numel(solid_face_keys)
        key = solid_face_keys{key_index};

        if face_count(key) == 1
            surface_faces{end + 1, 1} = ...
                face_oriented(key);
        end
    end

    surface_triangles = zeros(0, 3);

    for face_index = 1:numel(surface_faces)
        face_node_ids = surface_faces{face_index};

        face_rows = zeros(1, numel(face_node_ids));

        for j = 1:numel(face_node_ids)
            if ~isKey(id_to_row, face_node_ids(j))
                error('单元引用了不存在的节点 ID。');
            end

            face_rows(j) = id_to_row(face_node_ids(j));
        end

        if numel(face_rows) == 3
            surface_triangles(end + 1, :) = ...
                face_rows;

        elseif numel(face_rows) == 4
            surface_triangles(end + 1, :) = ...
                face_rows([1, 2, 3]);

            surface_triangles(end + 1, :) = ...
                face_rows([1, 3, 4]);
        end
    end

    if isempty(surface_triangles)
        error('没有生成有效的表面三角形。');
    end
end


function number_of_nodes = ...
    nodes_per_supported_element(element_type)

    if startsWith(element_type, 'S3')
        number_of_nodes = 3;

    elseif startsWith(element_type, 'S4')
        number_of_nodes = 4;

    elseif startsWith(element_type, 'C3D4')
        number_of_nodes = 4;

    elseif startsWith(element_type, 'C3D8')
        number_of_nodes = 8;

    else
        number_of_nodes = 0;
    end
end


function add_solid_faces( ...
    connectivity, local_faces, ...
    node_coords, id_to_row, ...
    face_count, face_oriented)
% 将实体单元各面加入计数器。
% 两个实体共享的面计数为 2，最后会被删除。

    element_rows = zeros(1, numel(connectivity));

    for j = 1:numel(connectivity)
        element_rows(j) = id_to_row(connectivity(j));
    end

    element_centroid = mean( ...
        node_coords(element_rows, :), 1);

    for face_index = 1:numel(local_faces)
        face_ids = connectivity(local_faces{face_index});

        face_rows = zeros(1, numel(face_ids));

        for j = 1:numel(face_ids)
            face_rows(j) = id_to_row(face_ids(j));
        end

        face_points = node_coords(face_rows, :);
        face_centroid = mean(face_points, 1);

        face_normal = cross( ...
            face_points(2, :) - face_points(1, :), ...
            face_points(3, :) - face_points(1, :));

        % 将实体外表面法向调整为从单元内部指向外部
        if dot( ...
                face_normal, ...
                face_centroid - element_centroid) < 0

            face_ids = face_ids([1, end:-1:2]);
        end

        sorted_ids = sort(face_ids);
        key = sprintf('%.0f_', sorted_ids);

        if isKey(face_count, key)
            face_count(key) = face_count(key) + 1;
        else
            face_count(key) = 1;
            face_oriented(key) = face_ids;
        end
    end
end
