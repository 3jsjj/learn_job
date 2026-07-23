% =========================================================================
% 曲面压力逆有限元重构（仅位移反演版）
%
% 核心逻辑：
%   1) 从 measured_filepath 读取测点 NodeID、U1、U2、U3；
%   2) 只使用测点位移作为观测量，反演目标曲面的压力；
%   3) CSV 第 5 列若存在，仅作为真实压力参考值，用于反演后比较，
%      绝不会作为已知载荷参与求解；
%   4) 默认先使用“均匀压力”模式验证整个反演链条，验证通过后可切换
%      到“节点压力”模式。
%
% 曲面压力到等效节点力：
%
%       f_i = pressure_sign * p_i * A_i * n_i
%
% 其中：
%   n_i：曲面节点单位法向；
%   A_i：节点分摊面积；
%   p_i：待反演压力；
%   pressure_sign = -1 表示正压力沿外法向反方向作用。
%
% 测点 CSV：
%   必需列：NodeID, U1, U2, U3
%   可选第 5 列：PressureReference，仅用于结果验证，不参与反演。
%
% 单位：
%   若刚度为 N/mm、坐标为 mm、位移为 mm，则压力为 N/mm^2，即 MPa。
%
% 本版本适用于当前壳单元刚度矩阵：
%   每节点 6 自由度，压力只装配到 U1/U2/U3 平移自由度。
% =========================================================================

clear;
clc;

%% 1. 用户设置

inp_filepath      = 'Zhimian_addnodes_cut.inp';
mtx_filepath      = 'Job-5-1_STIF2.mtx';
measured_filepath = 'Abaqus_Nodal_U_and_Pressure.csv';
fixed_filepath    = 'fixed_nodes.csv';

% 可选文件：每行 [X, Y, Z]，必须对应目标加载曲面的有限元节点。
% 若文件不存在，程序默认使用模型全部外表面节点。
virtual_coords_filepath = 'virtual_coords.csv';

% 当前 .mtx 为壳单元矩阵：每个节点 6 个自由度。
% 1-3 为平移 U1/U2/U3，4-6 为转动 UR1/UR2/UR3。
num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];

% 固定边界约束的自由度分量。
% 完全固支壳边界通常约束 1:6。
fixed_dof_components = 1:6;

% +1：正压力沿计算出的表面法向；
% -1：正压力沿表面法向反方向，外部压力压向结构时通常使用 -1。
pressure_sign = -1;

% 虚拟坐标到有限元节点的绝对匹配容差。
coordinate_tolerance = 1e-6;

% 反演模式：
%   'uniform'：整个目标曲面只有一个均匀压力未知量，建议首先使用；
%   'nodal'  ：每个虚拟节点一个压力未知量，使用 L-curve 正则化。
inversion_mode = 'uniform';

% 测点节点也允许作为压力未知节点。
% 因为测点位移是观测量，并不代表该节点压力已知。
exclude_measured_nodes_from_virtual = false;

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
    % 去掉了没有参与画三角形网格的nodes
    virtual_node_rows_default = unique(exterior_triangles(:));
    virtual_coords = node_coords(virtual_node_rows_default, :);

    warning(['未找到 virtual_coords.csv，', ...
             '当前使用模型全部外表面节点作为虚拟压力节点。']);
end

[virtual_node_rows, virtual_match_distance] = ...
    match_coordinates_to_nodes(virtual_coords, node_coords);
% 做一个筛选，节点和虚拟测点的位置不能隔太远
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


%% 4. 读取真实测点位移

measured_data = readmatrix(measured_filepath);
pressure_true = measured_data(:,5);
pressure_true = pressure_true(isfinite(pressure_true));

fprintf('\n========== 真实压力检查 ==========\n');
fprintf('压力数量：%d\n', numel(pressure_true));
fprintf('最小压力：%.9e\n', min(pressure_true));
fprintf('最大压力：%.9e\n', max(pressure_true));
fprintf('平均压力：%.9e\n', mean(pressure_true));
fprintf('标准差：%.9e\n', std(pressure_true));
fprintf('=================================\n');

% 反演只要求 NodeID、U1、U2、U3 四列。
if size(measured_data, 2) < 4
    error('测点文件至少需要四列：NodeID, U1, U2, U3。');
end

id_measured_raw = measured_data(:, 1);
U_measured_raw  = measured_data(:, 2:4);

% 如果恰好存在第 5 列，将其保存为真实压力参考值。
% 此列只用于反演完成后的比较，绝不参与载荷组装。
has_pressure_reference = size(measured_data, 2) >= 5 && ...
                         size(measured_data, 2) < 7;

pressure_reference_raw = nan(size(measured_data, 1), 1);

if has_pressure_reference
    pressure_reference_raw = measured_data(:, 5);
end

% 有效性筛选只依据 NodeID 和位移；参考压力允许为 NaN。
valid_measured_data = ...
    isfinite(id_measured_raw) & ...
    all(isfinite(U_measured_raw), 2);

id_measured_raw = ...
    id_measured_raw(valid_measured_data);

U_measured_raw = ...
    U_measured_raw(valid_measured_data, :);

pressure_reference_raw = ...
    pressure_reference_raw(valid_measured_data);

% 将测点 NodeID 映射到 INP 节点行号。
[measured_id_exists, measured_node_rows_raw] = ...
    ismember(id_measured_raw, node_ids);

if any(~measured_id_exists)
    warning('%d 个测点 NodeID 不在 inp 节点中，已剔除。', ...
        nnz(~measured_id_exists));
end

id_measured_raw = ...
    id_measured_raw(measured_id_exists);

U_measured_raw = ...
    U_measured_raw(measured_id_exists, :);

pressure_reference_raw = ...
    pressure_reference_raw(measured_id_exists);

measured_node_rows = ...
    measured_node_rows_raw(measured_id_exists);

if isempty(measured_node_rows)
    error('没有有效的真实测点位移。');
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
% 计算节点法向，以及每个节点所占的面积
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

% -------------------------------------------------------------------------
% 自动识别 MTX 节点编号与 INP 节点数据之间的关系
%
% node_ids/node_coords 的第 r 行描述同一个有限元节点。
% 当前 MTX 节点号可能是：
%   A. 直接使用 INP NodeID；
%   B. 使用节点在 node_ids/node_coords 中的行号。
%
% 构造：
%   node_row_to_mtx_id(r) = 第 r 个 INP 节点在 MTX 中的节点编号
% -------------------------------------------------------------------------
% 将node_i和node_j中的所有元素不重复的赋给mtx_node_ids
mtx_node_ids = unique([node_i; node_j]);
N_inp_nodes = numel(node_ids);

node_row_to_mtx_id = nan(N_inp_nodes, 1);

if all(ismember(mtx_node_ids, node_ids))
    % 情况 A：MTX 与 INP 直接使用相同 NodeID
    [~, mtx_node_rows] = ismember(mtx_node_ids, node_ids);
    node_row_to_mtx_id(mtx_node_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = INP NodeID';

elseif all(mtx_node_ids >= 1) && ...
       all(mtx_node_ids <= N_inp_nodes) && ...
       all(mtx_node_ids == round(mtx_node_ids))
    % 情况 B：MTX NodeID 等于 node_ids/node_coords 的行号
    mtx_node_rows = mtx_node_ids;
    node_row_to_mtx_id(mtx_node_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = node_ids/node_coords 行号';

else
    error([ ...
        '无法自动识别 MTX 节点编号与 INP 节点的对应关系。', newline, ...
        'MTX NodeID 既不是 INP NodeID，也不是有效节点行号。']);
end

fprintf('\n');
fprintf('============== 节点编号自动映射 ==============\n');
fprintf('映射模式：%s\n', node_mapping_mode);
fprintf('INP 节点数量：%d\n', N_inp_nodes);
fprintf('MTX 节点数量：%d\n', numel(mtx_node_ids));
fprintf('MTX NodeID 范围：[%d, %d]\n', ...
    min(mtx_node_ids), max(mtx_node_ids));
fprintf('INP 中有 MTX 映射的节点：%d\n', ...
    nnz(isfinite(node_row_to_mtx_id)));
fprintf('==============================================\n');

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
% nnz()是number of nonzeros
% triu是上半对角上的数据，K_INPUT是对象，1是对角线以上不包括对角线的意思，后面的tril是下半边，-1是不包括对角线的下半部分
has_upper = nnz(triu(K_input, 1)) > 0;
has_lower = nnz(tril(K_input, -1)) > 0;
% 两个条件有且仅有一个为真
if xor(has_upper, has_lower)
    %diag是提取主对角线的数据变成列，spdiags(元素，0=放在对角线，行数，列数)
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


%% 7. 由 fixed_nodes.csv 的 ID 自动取得坐标和 MTX 节点号

fixed_nodes_data = readmatrix(fixed_filepath);

fixed_nodes_input = unique( ...
    fixed_nodes_data(isfinite(fixed_nodes_data)));
fixed_nodes_input = fixed_nodes_input(:);

valid_fixed_id = ...
    fixed_nodes_input > 0 & ...
    abs(fixed_nodes_input - round(fixed_nodes_input)) < 1e-10;

if any(~valid_fixed_id)
    warning('%d 个固定节点条目不是正整数 NodeID，已剔除。', ...
        nnz(~valid_fixed_id));
end

fixed_nodes_input = round( ...
    fixed_nodes_input(valid_fixed_id));

% 第一步：固定 ID 在 node_ids 中定位
[fixed_found_in_inp, fixed_node_rows] = ...
    ismember(fixed_nodes_input, node_ids);

if any(~fixed_found_in_inp)
    missing_fixed_ids = ...
        fixed_nodes_input(~fixed_found_in_inp);

    print_count = min(10, numel(missing_fixed_ids));

    fprintf('无法在 node_ids 中找到的前 %d 个固定 ID：', ...
        print_count);
    fprintf('%d ', missing_fixed_ids(1:print_count));
    fprintf('\n');

    error([ ...
        'fixed_nodes.csv 中有 %d 个 ID 不在 INP 的 node_ids 中。', ...
        newline, ...
        'fixed_nodes.csv 必须使用与 inp 相同的节点标签。'], ...
        nnz(~fixed_found_in_inp));
end

% 第二步：直接从已读取的 node_coords 获取固定节点坐标
coords_fixed = node_coords(fixed_node_rows, :);

% 第三步：通过节点行号得到 MTX 使用的节点编号
fixed_nodes_mtx = ...
    node_row_to_mtx_id(fixed_node_rows);

fixed_has_mtx_mapping = isfinite(fixed_nodes_mtx);

fprintf('\n');
fprintf('============== 固定节点自动映射 ==============\n');
fprintf('fixed_nodes.csv 节点数量：%d\n', ...
    numel(fixed_nodes_input));
fprintf('在 INP node_ids 中找到：%d / %d\n', ...
    nnz(fixed_found_in_inp), numel(fixed_nodes_input));
fprintf('能映射到 MTX：%d / %d\n', ...
    nnz(fixed_has_mtx_mapping), numel(fixed_nodes_input));
fprintf('固定节点 X 范围：[%.9g, %.9g]\n', ...
    min(coords_fixed(:,1)), max(coords_fixed(:,1)));
fprintf('固定节点 Y 范围：[%.9g, %.9g]\n', ...
    min(coords_fixed(:,2)), max(coords_fixed(:,2)));
fprintf('固定节点 Z 范围：[%.9g, %.9g]\n', ...
    min(coords_fixed(:,3)), max(coords_fixed(:,3)));
fprintf('==============================================\n');

if any(~fixed_has_mtx_mapping)
    bad_fixed_ids = ...
        fixed_nodes_input(~fixed_has_mtx_mapping);

    print_count = min(10, numel(bad_fixed_ids));

    fprintf('没有 MTX 映射的前 %d 个固定 INP NodeID：', ...
        print_count);
    fprintf('%d ', bad_fixed_ids(1:print_count));
    fprintf('\n');

    error([ ...
        '部分固定节点存在于 INP，但没有出现在 MTX 中。', newline, ...
        '请确认 INP 和 MTX 来自同一个模型和同一次编号。']);
end

fixed_nodes_mtx = round(fixed_nodes_mtx);

% 壳单元完全固支：删除每节点 1～6 自由度
constrained_dofs_matrix = node_dof_numbers( ...
    fixed_nodes_mtx, ...
    num_dofs_per_node, ...
    fixed_dof_components);

constrained_dofs = constrained_dofs_matrix(:);

[nonzero_row, nonzero_col] = find(K_raw);
existing_dofs = unique([nonzero_row; nonzero_col]);

constrained_dofs_existing = intersect( ...
    constrained_dofs, existing_dofs);
% 减去了所有的固定约束的自由度标号的刚度值
active_dofs = setdiff( ...
    existing_dofs, constrained_dofs_existing);

K_global = K_raw(active_dofs, active_dofs);

fprintf('\n');
fprintf('刚度矩阵原始非零自由度：%d\n', ...
    numel(existing_dofs));
fprintf('固定节点数量：%d\n', ...
    numel(fixed_nodes_mtx));
fprintf('理论应删除的约束自由度：%d\n', ...
    numel(fixed_nodes_mtx) * numel(fixed_dof_components));
fprintf('实际删除的约束自由度：%d\n', ...
    numel(constrained_dofs_existing));
fprintf('活动自由度数量：%d\n', ...
    size(K_global, 1));

% 画出自动识别出的固定节点
figure( ...
    'Name', '固定边界自动识别检查', ...
    'Color', 'w', ...
    'Position', [100, 100, 900, 700]);

trisurf( ...
    exterior_triangles, ...
    node_coords(:,1), ...
    node_coords(:,2), ...
    node_coords(:,3), ...
    'FaceColor', [0.75, 0.80, 0.90], ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', [0.65, 0.65, 0.65]);

hold on;

scatter3( ...
    coords_fixed(:,1), ...
    coords_fixed(:,2), ...
    coords_fixed(:,3), ...
    35, 'r', 'filled');

axis equal;
grid on;
view(45, 30);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('红色节点应与 Abaqus 固定边缘一致');
legend('模型表面', '固定节点');
hold off;


%% 8. 将测点 INP NodeID 映射为 MTX NodeID，并过滤无效测点

measured_mtx_ids_raw = ...
    node_row_to_mtx_id(measured_node_rows);

measured_has_mtx_mapping = ...
    isfinite(measured_mtx_ids_raw);

measured_global_dofs = nan( ...
    numel(measured_mtx_ids_raw), ...
    numel(translation_dofs));

measured_global_dofs(measured_has_mtx_mapping, :) = ...
    node_dof_numbers( ...
        measured_mtx_ids_raw(measured_has_mtx_mapping), ...
        num_dofs_per_node, ...
        translation_dofs);

valid_measured_active = ...
    measured_has_mtx_mapping & ...
    all(ismember(measured_global_dofs, active_dofs), 2);

if any(~valid_measured_active)
    warning('%d 个测点没有 MTX 映射、属于固定边界或无活动自由度，已剔除。', ...
        nnz(~valid_measured_active));
end

id_measured = ...
    id_measured_raw(valid_measured_active);

measured_node_rows = ...
    measured_node_rows(valid_measured_active);

measured_mtx_ids = ...
    measured_mtx_ids_raw(valid_measured_active);

U_measured_mat = ...
    U_measured_raw(valid_measured_active, :);

pressure_reference = ...
    pressure_reference_raw(valid_measured_active);

measured_global_dofs = ...
    measured_global_dofs(valid_measured_active, :);

[~, dof_measured_matrix] = ismember( ...
    measured_global_dofs, active_dofs);

dof_measured = reshape( ...
    dof_measured_matrix', [], 1);

U_measured = reshape( ...
    U_measured_mat', [], 1);

if isempty(dof_measured)
    error('过滤后没有有效测点自由度。');
end

%% 9. 反演时不使用真实压力作为已知载荷

N_total_dof = size(K_global, 1);

% 这里只使用测点位移反演压力。
% 即使 measured CSV 第 5 列存在真实压力，也不把它加入 F_known。
F_known = sparse(N_total_dof, 1);

fprintf('\n');
fprintf('反演模式：只使用测点位移，已知载荷向量 F_known = 0。\n');

%% 10. 过滤虚拟压力节点并构造压力到节点力的映射 G

% 输出时保留 INP NodeID；自由度计算使用 MTX NodeID
id_virtual = node_ids(virtual_node_rows);
virtual_mtx_ids = node_row_to_mtx_id(virtual_node_rows);

virtual_has_mtx_mapping = isfinite(virtual_mtx_ids);

virtual_global_dofs = nan( ...
    numel(id_virtual), numel(translation_dofs));

virtual_global_dofs(virtual_has_mtx_mapping, :) = ...
    node_dof_numbers( ...
        virtual_mtx_ids(virtual_has_mtx_mapping), ...
        num_dofs_per_node, ...
        translation_dofs);

% 将各种过滤原因分别统计
virtual_has_stiffness = ...
    virtual_has_mtx_mapping & ...
    all(ismember(virtual_global_dofs, existing_dofs), 2);

virtual_is_fixed = ismember( ...
    virtual_node_rows, fixed_node_rows);

virtual_is_active = ...
    virtual_has_mtx_mapping & ...
    all(ismember(virtual_global_dofs, active_dofs), 2);

virtual_is_measured = ismember( ...
    virtual_node_rows, measured_node_rows);

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
virtual_mtx_ids   = virtual_mtx_ids(valid_virtual_active);
virtual_node_rows = virtual_node_rows(valid_virtual_active);
virtual_coords    = virtual_coords(valid_virtual_active, :);
virtual_global_dofs = virtual_global_dofs(valid_virtual_active, :);

[~, dof_virtual_matrix] = ismember( ...
    virtual_global_dofs, active_dofs);

N_virtual_nodes = numel(id_virtual);

%% 检查真实压力节点与反演曲面的对应关系

pressure_reference_ids = measured_data(:,1);
pressure_reference_values = measured_data(:,5);

valid_pressure_reference = ...
    isfinite(pressure_reference_ids) & ...
    isfinite(pressure_reference_values);

pressure_reference_ids = ...
    pressure_reference_ids(valid_pressure_reference);

pressure_reference_values = ...
    pressure_reference_values(valid_pressure_reference);

% 真实压力节点是否属于反演节点
[pressure_is_virtual, pressure_virtual_index] = ...
    ismember(pressure_reference_ids, id_virtual);

% 反演节点中哪些具有真实压力
[virtual_has_true_pressure, virtual_true_index] = ...
    ismember(id_virtual, pressure_reference_ids);

fprintf('\n');
fprintf('========== 压力节点对应检查 ==========\n');

fprintf('真实压力数据数量：%d\n', ...
    numel(pressure_reference_ids));

fprintf('反演压力节点数量：%d\n', ...
    N_virtual_nodes);

fprintf('真实压力点属于反演节点：%d / %d\n', ...
    nnz(pressure_is_virtual), ...
    numel(pressure_reference_ids));

fprintf('反演节点具有真实压力：%d / %d\n', ...
    nnz(virtual_has_true_pressure), ...
    N_virtual_nodes);

fprintf('覆盖比例：%.2f %%\n', ...
    100 * nnz(virtual_has_true_pressure) / ...
    max(N_virtual_nodes, 1));

if any(~pressure_is_virtual)
    fprintf('没有进入反演曲面的真实压力 NodeID：\n');
    disp(pressure_reference_ids(~pressure_is_virtual));
end

fprintf('======================================\n');



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

% 说明一下G (A,B) VALUE=force_per_unit_pressure, 这里的A是在自由度编号，B是节点编号
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


%% 12. 将测点位移直接作为反演目标

% 本版本没有已知压力载荷，因此不需要扣除基准位移。
U_from_known_force = zeros(N_total_dof, 1);
U_measured_baseline = zeros(size(U_measured));

Delta_U = U_measured;

fprintf('\n');
fprintf('输入测点位移范数：%.6e\n', norm(U_measured));
fprintf('反演目标位移范数：%.6e\n', norm(Delta_U));

%% 13. 使用测点位移反演曲面压力

fprintf('\n');
fprintf('当前压力反演模式：%s\n', inversion_mode);

switch lower(inversion_mode)

    case 'uniform'
        % ---------------------------------------------------------------
        % 均匀压力模式：整个目标曲面只有一个压力未知量 p_uniform。
        % pressure_virtual = p_uniform * ones(N_virtual_nodes, 1)
        % 推荐先用此模式验证法向、面积、单位、边界和刚度矩阵。
        % ---------------------------------------------------------------
        uniform_pressure_basis = ...
            ones(N_virtual_nodes, 1);

        H_uniform = ...
            H * uniform_pressure_basis;

        uniform_denominator = ...
            H_uniform' * H_uniform;

        if uniform_denominator <= eps
            error([ ...
                '均匀压力对测点位移几乎没有灵敏度。', newline, ...
                '请检查加载曲面、测点位置、边界条件和压力方向。']);
        end

        % 一个未知量的最小二乘解。
        p_uniform = ...
            (H_uniform' * Delta_U) / uniform_denominator;

        pressure_virtual = ...
            p_uniform * uniform_pressure_basis;

        lambda_opt = 0;

        fprintf('均匀压力反演结果：%.9e\n', p_uniform);

    case 'nodal'
        % ---------------------------------------------------------------
        % 节点压力模式：每个虚拟节点一个压力未知量。
        % 该问题通常欠定或病态，使用 L-curve 选择 Tikhonov 参数。
        % ---------------------------------------------------------------
        disp('正在使用 L-curve 法寻找节点压力正则化参数...');

        A = full(H * H');
        A = (A + A') / 2;

        [Q, D] = eig(A);
        d = max(real(diag(D)), 0);
        d_max = max(d);

        if isempty(d) || d_max <= eps
            error('灵敏度矩阵 H 几乎为零，无法进行压力反演。');
        end

        u_projected = Q' * Delta_U;

        n_lambdas = 100;
        lambda_vec = logspace(-12, 0, n_lambdas) * d_max;
        t_vec = log10(lambda_vec);

        eta = zeros(1, n_lambdas);
        rho = zeros(1, n_lambdas);

        for i = 1:n_lambdas
            lam = lambda_vec(i);
            w = u_projected ./ (d + lam);

            eta(i) = sqrt(sum(d .* abs(w).^2));
            rho(i) = sqrt(sum(abs( ...
                lam .* u_projected ./ (d + lam)).^2));
        end

        rho_safe = max(rho, eps);
        eta_safe = max(eta, eps);

        x = log(rho_safe);
        y = log(eta_safe);

        dx_dt = gradient(x, t_vec);
        dy_dt = gradient(y, t_vec);
        d2x_dt2 = gradient(dx_dt, t_vec);
        d2y_dt2 = gradient(dy_dt, t_vec);

        curvature = ...
            (dx_dt .* d2y_dt2 - dy_dt .* d2x_dt2) ./ ...
            max((dx_dt.^2 + dy_dt.^2).^(1.5), eps);

        curvature_score = abs(curvature);
        curvature_score(~isfinite(curvature_score)) = -inf;

        valid_idx = 5:(n_lambdas - 5);
        [~, max_idx_offset] = ...
            max(curvature_score(valid_idx));

        opt_idx = valid_idx(max_idx_offset);
        lambda_opt = lambda_vec(opt_idx);

        fprintf('L-curve 寻优完成，最优 lambda = %.6e\n', ...
            lambda_opt);

        figure( ...
            'Name', 'L-Curve Analysis', ...
            'Color', 'w', ...
            'Position', [150, 150, 1000, 400]);

        subplot(1, 2, 1);
        loglog(rho_safe, eta_safe, 'b-', 'LineWidth', 2);
        hold on;
        loglog( ...
            rho_safe(opt_idx), eta_safe(opt_idx), ...
            'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        grid on;
        xlabel('残差范数 ||Hp - U_m||');
        ylabel('解的范数 ||p||');
        title('L-Curve');
        legend( ...
            'L 曲线', ...
            sprintf('选取点 (\\lambda = %.2e)', lambda_opt), ...
            'Location', 'northeast');

        subplot(1, 2, 2);
        semilogx(lambda_vec, curvature_score, ...
            'k-', 'LineWidth', 1.5);
        hold on;
        semilogx( ...
            lambda_opt, curvature_score(opt_idx), ...
            'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        grid on;
        xlabel('正则化参数 \lambda');
        ylabel('|曲率|');
        title('L-Curve 曲率');

        dual_matrix = ...
            A + lambda_opt * speye(N_measured_dof);

        pressure_virtual = ...
            H' * (dual_matrix \ Delta_U);

    otherwise
        error([ ...
            '未知 inversion_mode：', inversion_mode, newline, ...
            '请设置为 uniform 或 nodal。']);
end

%% 14. 根据反演压力计算完整位移和虚拟点位移

F_pressure_global = ...
    G * pressure_virtual;

% 本版本没有其他已知载荷。
F_total_global = ...
    F_pressure_global;

U_reconstructed_global = ...
    K_global \ F_total_global;

dof_virtual = reshape( ...
    dof_virtual_matrix', [], 1);

U_virtual_vector = ...
    U_reconstructed_global(dof_virtual);

U_virtual_3D = reshape( ...
    U_virtual_vector, 3, [])';

% 压力产生的三维等效节点力，仅用于检查。
F_virtual_equivalent = bsxfun( ...
    @times, force_per_unit_pressure, pressure_virtual);

%% 15. 误差、可辨识性和参考压力诊断

U_measured_predicted = ...
    U_reconstructed_global(dof_measured);

U_measured_predicted_mat = reshape( ...
    U_measured_predicted, 3, [])';

U_measured_error_mat = ...
    U_measured_predicted_mat - U_measured_mat;

relative_displacement_error = ...
    norm(U_measured_predicted - U_measured) / ...
    max(norm(U_measured), eps);

equilibrium_residual = ...
    norm(K_global * U_reconstructed_global - F_total_global) / ...
    max(norm(F_total_global), eps);

% 检查灵敏度矩阵的可辨识性。
singular_values = svd(full(H), 'econ');

if isempty(singular_values)
    effective_rank = 0;
else
    rank_tolerance = ...
        max(size(H)) * eps(max(singular_values));

    effective_rank = ...
        nnz(singular_values > rank_tolerance);
end

fprintf('\n');
fprintf('================ 重构诊断 ================\n');
fprintf('有效测点数量：%d\n', numel(id_measured));
fprintf('测量自由度数量：%d\n', N_measured_dof);
fprintf('虚拟压力节点数量：%d\n', N_virtual_nodes);
fprintf('灵敏度矩阵有效秩：%d\n', effective_rank);
fprintf('测点位移相对拟合误差：%.6e\n', ...
    relative_displacement_error);
fprintf('全局平衡相对残差：%.6e\n', ...
    equilibrium_residual);
fprintf('最小重构压力：%.6e\n', min(pressure_virtual));
fprintf('最大重构压力：%.6e\n', max(pressure_virtual));
fprintf('==========================================\n');

% 输出每个测点的输入位移、预测位移和误差。
measured_displacement_table = table( ...
    id_measured(:), ...
    U_measured_mat(:, 1), ...
    U_measured_predicted_mat(:, 1), ...
    U_measured_error_mat(:, 1), ...
    U_measured_mat(:, 2), ...
    U_measured_predicted_mat(:, 2), ...
    U_measured_error_mat(:, 2), ...
    U_measured_mat(:, 3), ...
    U_measured_predicted_mat(:, 3), ...
    U_measured_error_mat(:, 3), ...
    'VariableNames', { ...
        'NodeID', ...
        'U1_Measured', 'U1_Predicted', 'U1_Error', ...
        'U2_Measured', 'U2_Predicted', 'U2_Error', ...
        'U3_Measured', 'U3_Predicted', 'U3_Error'});

writetable( ...
    measured_displacement_table, ...
    'Measured_vs_Predicted_Displacement.csv');

% 将可选的真实压力参考值映射到反演节点，只用于比较。
pressure_reference_on_virtual = ...
    nan(N_virtual_nodes, 1);

[reference_node_exists, reference_measured_index] = ...
    ismember(id_virtual, id_measured);

reference_candidate = ...
    reference_node_exists;

pressure_reference_on_virtual(reference_candidate) = ...
    pressure_reference( ...
        reference_measured_index(reference_candidate));

valid_pressure_reference = ...
    isfinite(pressure_reference_on_virtual);

pressure_error_on_virtual = ...
    nan(N_virtual_nodes, 1);

pressure_relative_error_on_virtual = ...
    nan(N_virtual_nodes, 1);

if any(valid_pressure_reference)
    pressure_error_on_virtual(valid_pressure_reference) = ...
        pressure_virtual(valid_pressure_reference) - ...
        pressure_reference_on_virtual(valid_pressure_reference);

    pressure_relative_error_on_virtual(valid_pressure_reference) = ...
        abs(pressure_error_on_virtual(valid_pressure_reference)) ./ ...
        max(abs(pressure_reference_on_virtual(valid_pressure_reference)), eps);

    pressure_reference_rmse = sqrt(mean( ...
        pressure_error_on_virtual(valid_pressure_reference).^2));

    pressure_reference_relative_L2 = ...
        norm(pressure_error_on_virtual(valid_pressure_reference)) / ...
        max(norm(pressure_reference_on_virtual(valid_pressure_reference)), eps);

    fprintf('\n');
    fprintf('参考压力匹配节点数量：%d\n', ...
        nnz(valid_pressure_reference));
    fprintf('参考压力 RMSE：%.6e\n', ...
        pressure_reference_rmse);
    fprintf('参考压力相对 L2 误差：%.2f %%\n', ...
        100 * pressure_reference_relative_L2);
    fprintf(['注意：参考压力只用于这里的结果比较，', ...
             '没有参与反演求解。\n']);

    pressure_reference_comparison_table = table( ...
        id_virtual(valid_pressure_reference), ...
        pressure_reference_on_virtual(valid_pressure_reference), ...
        pressure_virtual(valid_pressure_reference), ...
        pressure_error_on_virtual(valid_pressure_reference), ...
        pressure_relative_error_on_virtual(valid_pressure_reference), ...
        'VariableNames', { ...
            'NodeID', ...
            'Pressure_Reference', ...
            'Pressure_Reconstructed', ...
            'Pressure_Error', ...
            'Pressure_RelativeError'});

    writetable( ...
        pressure_reference_comparison_table, ...
        'Pressure_Reference_vs_Reconstructed.csv');
end

%% 16. 输出 CSV

result_table = table( ...
    id_virtual(:), ...
    virtual_mtx_ids(:), ...
    virtual_coords(:, 1), ...
    virtual_coords(:, 2), ...
    virtual_coords(:, 3), ...
    virtual_normals(:, 1), ...
    virtual_normals(:, 2), ...
    virtual_normals(:, 3), ...
    virtual_areas(:), ...
    pressure_virtual(:), ...
    pressure_reference_on_virtual(:), ...
    pressure_error_on_virtual(:), ...
    pressure_relative_error_on_virtual(:), ...
    U_virtual_3D(:, 1), ...
    U_virtual_3D(:, 2), ...
    U_virtual_3D(:, 3), ...
    F_virtual_equivalent(:, 1), ...
    F_virtual_equivalent(:, 2), ...
    F_virtual_equivalent(:, 3), ...
    'VariableNames', { ...
        'INP_NodeID', ...
        'MTX_NodeID', ...
        'X', 'Y', 'Z', ...
        'NormalX', 'NormalY', 'NormalZ', ...
        'NodalArea', ...
        'Pressure_Reconstructed', ...
        'Pressure_Reference', ...
        'Pressure_Error', ...
        'Pressure_RelativeError', ...
        'U1_Reconstructed', ...
        'U2_Reconstructed', ...
        'U3_Reconstructed', ...
        'EquivalentForceX', ...
        'EquivalentForceY', ...
        'EquivalentForceZ'});

writetable( ...
    result_table, ...
    'Displacement_Based_Pressure_Reconstruction.csv');

disp('结果已保存：Displacement_Based_Pressure_Reconstruction.csv');

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
title(['仅位移反演曲面压力：', inversion_mode]);

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

%% 17. 绘制虚拟节点 U1、U2、U3 位移云图

% 位移单位，请根据模型单位修改
displacement_unit = 'mm';

% exterior_triangles 中保存的是 node_coords 的行号。
% 如果你的变量名是 surface_triangles，请改为：
% exterior_triangles = surface_triangles;

% 基本尺寸检查
if size(U_virtual_3D, 1) ~= numel(virtual_node_rows)
    error(['U_virtual_3D 的行数与 virtual_node_rows 数量不一致。', ...
        '无法建立节点位移和节点位置的对应关系。']);
end

if size(U_virtual_3D, 2) < 3
    error('U_virtual_3D 至少需要包含 U1、U2、U3 三列。');
end

%% 将外表面三角形的全模型节点行号，映射成虚拟节点局部行号

% 第 r 个全模型节点对应虚拟节点数组中的第几个节点。
% 0 表示该节点不在当前虚拟节点集合中。
node_row_to_virtual_row = zeros(size(node_coords, 1), 1);

node_row_to_virtual_row(virtual_node_rows) = ...
    (1:numel(virtual_node_rows))';

% 将三角形中的 node_coords 行号转换为 virtual_coords 行号
virtual_triangles = ...
    node_row_to_virtual_row(exterior_triangles);

% 只有三个顶点都属于虚拟节点集合的三角形才能绘制
valid_virtual_triangles = all(virtual_triangles > 0, 2);

virtual_triangles = ...
    virtual_triangles(valid_virtual_triangles, :);

fprintf('\n位移云图使用的有效三角形数量：%d / %d\n', ...
    size(virtual_triangles, 1), ...
    size(exterior_triangles, 1));

if isempty(virtual_triangles)
    warning(['没有完整保留的虚拟节点三角形。', ...
        '程序将自动使用散点云图。']);
end

%% 三个位移分量分别绘图

component_names = {'U1', 'U2', 'U3'};
component_titles = { ...
    '重构位移云图 U1', ...
    '重构位移云图 U2', ...
    '重构位移云图 U3'};

for component_index = 1:3

    displacement_value = ...
        U_virtual_3D(:, component_index);

    figure( ...
        'Name', component_titles{component_index}, ...
        'Color', 'w', ...
        'Position', [100, 100, 850, 650]);

    if ~isempty(virtual_triangles)

        % 基于外表面三角形绘制连续云图
        trisurf( ...
            virtual_triangles, ...
            virtual_coords(:, 1), ...
            virtual_coords(:, 2), ...
            virtual_coords(:, 3), ...
            displacement_value, ...
            'FaceColor', 'interp', ...
            'EdgeColor', 'none');

    else

        % 没有完整三角形时，退化为节点散点云图
        scatter3( ...
            virtual_coords(:, 1), ...
            virtual_coords(:, 2), ...
            virtual_coords(:, 3), ...
            40, ...
            displacement_value, ...
            'filled', ...
            'MarkerEdgeColor', 'none');
    end

    axis equal;
    axis tight;
    grid on;
    box on;
    view(45, 30);

    xlabel('X 坐标');
    ylabel('Y 坐标');
    zlabel('Z 坐标');

    title(component_titles{component_index}, ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    colormap(jet(256));

    cb = colorbar;
    ylabel(cb, ...
        sprintf('%s 位移 (%s)', ...
        component_names{component_index}, ...
        displacement_unit), ...
        'FontSize', 11);

    % 令颜色范围关于 0 对称，便于观察正负位移
    finite_values = displacement_value( ...
        isfinite(displacement_value));

    if ~isempty(finite_values)
        color_limit = max(abs(finite_values));

        if color_limit > 0
            caxis([-color_limit, color_limit]);
        end
    end

    set(gca, ...
        'FontSize', 11, ...
        'DataAspectRatio', [1, 1, 1]);

    % 可选：保存 300 dpi 图片
    output_image_name = sprintf( ...
        'Reconstructed_%s_Cloud.png', ...
        component_names{component_index});

    exportgraphics( ...
        gcf, ...
        output_image_name, ...
        'Resolution', 300);

    fprintf('已保存：%s\n', output_image_name);
end


%% 18. 重新计算完整外表面位移，用于和 Abaqus 对比

% 反演后的总活动载荷
F_total_reconstructed = ...
    G * pressure_virtual;

% 反演载荷产生的全部活动自由度位移
U_full_reconstructed = ...
    K_global \ F_total_reconstructed;

% 获取全部外表面节点行号
surface_node_rows = unique(exterior_triangles(:));

surface_node_ids = node_ids(surface_node_rows);
surface_mtx_ids = node_row_to_mtx_id(surface_node_rows);

N_surface_nodes = numel(surface_node_rows);

% 初始化外表面位移
% 固定节点默认位移为 0；
% 无法映射的节点后面设置为 NaN
U_surface_3D = zeros(N_surface_nodes, 3);

surface_mapping_valid = isfinite(surface_mtx_ids);

surface_global_dofs = nan(N_surface_nodes, 3);

surface_global_dofs(surface_mapping_valid, :) = ...
    node_dof_numbers( ...
        surface_mtx_ids(surface_mapping_valid), ...
        num_dofs_per_node, ...
        translation_dofs);

% 检查三个平移自由度是否都属于活动自由度
surface_active = ...
    surface_mapping_valid & ...
    all(ismember(surface_global_dofs, active_dofs), 2);

% 将绝对自由度编号转换为 K_global 局部编号
[~, surface_local_dofs] = ...
    ismember(surface_global_dofs(surface_active, :), ...
             active_dofs);

% 提取三个方向位移
U_surface_3D(surface_active, :) = ...
    reshape( ...
        U_full_reconstructed( ...
            reshape(surface_local_dofs', [], 1)), ...
        3, [])';

% 无法进行节点映射的点设为 NaN
U_surface_3D(~surface_mapping_valid, :) = NaN;

% 外表面节点坐标
surface_coords = node_coords(surface_node_rows, 1:3);

%% 将原始外表面三角形转换为 surface_coords 的局部编号

node_row_to_surface_row = zeros(size(node_coords, 1), 1);

node_row_to_surface_row(surface_node_rows) = ...
    (1:N_surface_nodes)';

surface_triangles_local = ...
    node_row_to_surface_row(exterior_triangles);

valid_surface_triangles = ...
    all(surface_triangles_local > 0, 2);

surface_triangles_local = ...
    surface_triangles_local(valid_surface_triangles, :);

%% 分别绘制完整外表面的 U1、U2、U3

component_names = {'U1', 'U2', 'U3'};

for component_index = 1:3

    displacement_value = ...
        U_surface_3D(:, component_index);

    figure( ...
        'Name', ['完整表面 ', component_names{component_index}], ...
        'Color', 'w', ...
        'Position', [100, 100, 850, 650]);

    trisurf( ...
        surface_triangles_local, ...
        surface_coords(:, 1), ...
        surface_coords(:, 2), ...
        surface_coords(:, 3), ...
        displacement_value, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none');

    axis equal;
    axis tight;
    grid on;
    box on;
    view(45, 30);

    xlabel('X 坐标');
    ylabel('Y 坐标');
    zlabel('Z 坐标');

    title( ...
        sprintf('完整外表面重构位移 %s', ...
        component_names{component_index}), ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    colormap(jet(256));

    cb = colorbar;
    ylabel(cb, ...
        sprintf('%s 位移 (%s)', ...
        component_names{component_index}, ...
        displacement_unit));

    finite_values = displacement_value( ...
        isfinite(displacement_value));

    if ~isempty(finite_values)
        color_limit = max(abs(finite_values));

        if color_limit > 0
            caxis([-color_limit, color_limit]);
        end
    end

    output_image_name = sprintf( ...
        'Full_Surface_Reconstructed_%s.png', ...
        component_names{component_index});

    exportgraphics( ...
        gcf, output_image_name, ...
        'Resolution', 300);
end

fprintf('\n========== 均匀压力数量级诊断 ==========\n');

fprintf('测点输入位移最小值：%.6e\n', ...
    min(U_measured));

fprintf('测点输入位移最大值：%.6e\n', ...
    max(U_measured));

fprintf('测点输入位移范数：%.6e\n', ...
    norm(U_measured));

fprintf('单位均匀压力产生的测点位移最小值：%.6e\n', ...
    min(H_uniform));

fprintf('单位均匀压力产生的测点位移最大值：%.6e\n', ...
    max(H_uniform));

fprintf('单位均匀压力产生的测点位移范数：%.6e\n', ...
    norm(H_uniform));

fprintf('位移投影分子 H_uniform''*U：%.6e\n', ...
    H_uniform' * U_measured);

fprintf('灵敏度平方 H_uniform''*H_uniform：%.6e\n', ...
    H_uniform' * H_uniform);

fprintf('反演均匀压力：%.6e\n', ...
    p_uniform);

fprintf('=========================================\n');

U_predicted_uniform = ...
    H_uniform * p_uniform;

uniform_relative_error = ...
    norm(U_predicted_uniform - U_measured) / ...
    max(norm(U_measured), eps);

uniform_cosine_similarity = ...
    dot(U_predicted_uniform, U_measured) / ...
    max(norm(U_predicted_uniform) * ...
        norm(U_measured), eps);

fprintf('均匀压力预测位移范数：%.6e\n', ...
    norm(U_predicted_uniform));

fprintf('均匀压力位移相对误差：%.2f %%\n', ...
    100 * uniform_relative_error);

fprintf('预测与测量位移余弦相似度：%.6f\n', ...
    uniform_cosine_similarity);
%% 检查测量位移与单位压力响应的分量对应关系

U_measured_3D = reshape( ...
    U_measured, 3, [])';

H_uniform_3D = reshape( ...
    H_uniform, 3, [])';

%% 枚举检查轴交换和方向反转

axis_permutations = perms(1:3);

best_cosine = -inf;
best_error = inf;
best_permutation = [];
best_signs = [];
best_pressure = NaN;
best_U_transformed = [];

for permutation_index = 1:size(axis_permutations,1)

    current_permutation = ...
        axis_permutations(permutation_index,:);

    for sign_1 = [-1, 1]
        for sign_2 = [-1, 1]
            for sign_3 = [-1, 1]

                current_signs = ...
                    [sign_1, sign_2, sign_3];

                % 将测量位移重新排列并改变方向
                U_test = ...
                    U_measured_3D(:,current_permutation);

                U_test = ...
                    U_test .* current_signs;

                % 当前坐标变换下的最佳均匀压力比例
                p_test = ...
                    sum(H_uniform_3D(:) .* U_test(:)) / ...
                    max(sum(H_uniform_3D(:).^2), eps);

                U_test_predicted = ...
                    p_test * H_uniform_3D;

                cosine_test = ...
                    sum(U_test_predicted(:) .* U_test(:)) / ...
                    max(norm(U_test_predicted(:)) * ...
                        norm(U_test(:)), eps);

                error_test = ...
                    norm(U_test_predicted - U_test, 'fro') / ...
                    max(norm(U_test, 'fro'), eps);

                if cosine_test > best_cosine

                    best_cosine = cosine_test;
                    best_error = error_test;
                    best_permutation = current_permutation;
                    best_signs = current_signs;
                    best_pressure = p_test;
                    best_U_transformed = U_test;

                end
            end
        end
    end
end

fprintf('\n');
fprintf('========== 轴交换与方向反转检查 ==========\n');

fprintf('最佳列排列：[%d %d %d]\n', ...
    best_permutation(1), ...
    best_permutation(2), ...
    best_permutation(3));

fprintf('最佳方向符号：[%+d %+d %+d]\n', ...
    best_signs(1), ...
    best_signs(2), ...
    best_signs(3));

fprintf('最佳拟合压力：%.6e\n', ...
    best_pressure);

fprintf('变换后余弦相似度：%.6f\n', ...
    best_cosine);

fprintf('变换后相对误差：%.2f %%\n', ...
    100 * best_error);

fprintf('==========================================\n');





component_correlation = zeros(3, 3);
component_best_scale  = zeros(3, 3);

for measured_component = 1:3

    u_component = ...
        U_measured_3D(:, measured_component);

    for predicted_component = 1:3

        h_component = ...
            H_uniform_3D(:, predicted_component);

        % 余弦相似度
        component_correlation( ...
            measured_component, ...
            predicted_component) = ...
            dot(u_component, h_component) / ...
            max(norm(u_component) * ...
                norm(h_component), eps);

        % 最佳比例：
        % u_component ≈ scale * h_component
        component_best_scale( ...
            measured_component, ...
            predicted_component) = ...
            dot(h_component, u_component) / ...
            max(dot(h_component, h_component), eps);
    end
end

fprintf('\n');
fprintf('========== 位移分量对应关系检查 ==========\n');
fprintf('相关系数矩阵：\n');
fprintf('行 = 测量 U1/U2/U3，列 = 预测 U1/U2/U3\n');
disp(component_correlation);

fprintf('各分量最佳拟合压力系数：\n');
fprintf('行 = 测量 U1/U2/U3，列 = 预测 U1/U2/U3\n');
disp(component_best_scale);

fprintf('测量位移各分量范数：\n');
disp(vecnorm(U_measured_3D, 2, 1));

fprintf('单位压力响应各分量范数：\n');
disp(vecnorm(H_uniform_3D, 2, 1));

fprintf('==========================================\n');

%% 检查测量位移与预测位移是否存在整体坐标系旋转

% U_measured_3D：测量位移，尺寸 N×3
% H_uniform_3D：单位均匀压力产生的位移，尺寸 N×3

cross_matrix = ...
    H_uniform_3D' * U_measured_3D;

[rotation_U, ~, rotation_V] = ...
    svd(cross_matrix);

% 构造不含镜像的正交旋转矩阵
rotation_correction = eye(3);
rotation_correction(3,3) = ...
    sign(det(rotation_U * rotation_V'));

rotation_matrix = ...
    rotation_U * rotation_correction * rotation_V';

% 把预测位移旋转到测量坐标系
H_uniform_rotated = ...
    H_uniform_3D * rotation_matrix;

% 旋转后最优统一压力
p_after_rotation = ...
    sum(H_uniform_rotated(:) .* U_measured_3D(:)) / ...
    max(sum(H_uniform_rotated(:).^2), eps);

% 旋转并缩放后的预测位移
U_predicted_after_rotation = ...
    p_after_rotation * H_uniform_rotated;

rotation_relative_error = ...
    norm(U_predicted_after_rotation - U_measured_3D, 'fro') / ...
    max(norm(U_measured_3D, 'fro'), eps);

rotation_cosine_similarity = ...
    sum(U_predicted_after_rotation(:) .* U_measured_3D(:)) / ...
    max(norm(U_predicted_after_rotation(:)) * ...
        norm(U_measured_3D(:)), eps);

fprintf('\n');
fprintf('========== 整体坐标旋转检查 ==========\n');

fprintf('最优旋转矩阵：\n');
disp(rotation_matrix);

fprintf('旋转矩阵行列式：%.6f\n', ...
    det(rotation_matrix));

fprintf('旋转后的最优均匀压力：%.6e\n', ...
    p_after_rotation);

fprintf('旋转后的预测位移范数：%.6e\n', ...
    norm(U_predicted_after_rotation, 'fro'));

fprintf('旋转后的位移相对误差：%.2f %%\n', ...
    100 * rotation_relative_error);

fprintf('旋转后的余弦相似度：%.6f\n', ...
    rotation_cosine_similarity);

fprintf('======================================\n');

%% 检查节点位移模长是否一致
% 位移模长不受坐标系旋转影响

U_measured_magnitude = ...
    vecnorm(U_measured_3D, 2, 2);

H_uniform_magnitude = ...
    vecnorm(H_uniform_3D, 2, 2);

% 模长之间的相关系数
magnitude_corr_matrix = ...
    corrcoef(U_measured_magnitude, H_uniform_magnitude);

if size(magnitude_corr_matrix, 1) == 2
    magnitude_correlation = ...
        magnitude_corr_matrix(1,2);
else
    magnitude_correlation = NaN;
end

% 使用模长单独拟合一个压力比例
p_magnitude = ...
    dot(H_uniform_magnitude, U_measured_magnitude) / ...
    max(dot(H_uniform_magnitude, H_uniform_magnitude), eps);

U_magnitude_predicted = ...
    p_magnitude * H_uniform_magnitude;

magnitude_relative_error = ...
    norm(U_magnitude_predicted - U_measured_magnitude) / ...
    max(norm(U_measured_magnitude), eps);

fprintf('\n');
fprintf('========== 位移模长检查 ==========\n');

fprintf('真实测点位移模长范围：[%.6e, %.6e]\n', ...
    min(U_measured_magnitude), ...
    max(U_measured_magnitude));

fprintf('单位压力响应模长范围：[%.6e, %.6e]\n', ...
    min(H_uniform_magnitude), ...
    max(H_uniform_magnitude));

fprintf('位移模长相关系数：%.6f\n', ...
    magnitude_correlation);

fprintf('基于位移模长拟合的压力：%.6e\n', ...
    p_magnitude);

fprintf('位移模长相对拟合误差：%.2f %%\n', ...
    100 * magnitude_relative_error);

fprintf('==================================\n');

magnitude_table = table( ...
    id_measured(:), ...
    U_measured_magnitude, ...
    H_uniform_magnitude, ...
    U_magnitude_predicted, ...
    'VariableNames', { ...
        'NodeID', ...
        'MeasuredDisplacementMagnitude', ...
        'UnitPressureDisplacementMagnitude', ...
        'FittedDisplacementMagnitude'});

writetable( ...
    magnitude_table, ...
    'Displacement_Magnitude_Check.csv');

%% 压力单位数量级检查

pressure_rms_raw = ...
    sqrt(mean(pressure_true.^2));

fprintf('\n');
fprintf('========== 压力单位数量级检查 ==========\n');

fprintf('原始压力平均值：%.6e\n', ...
    mean(pressure_true));

fprintf('原始压力绝对值平均值：%.6e\n', ...
    mean(abs(pressure_true)));

fprintf('原始压力 RMS：%.6e\n', ...
    pressure_rms_raw);

fprintf('\n假设原始数据单位为 kPa：\n');
fprintf('换算为 MPa 后平均值：%.6e\n', ...
    mean(pressure_true) / 1000);

fprintf('换算为 MPa 后绝对值平均值：%.6e\n', ...
    mean(abs(pressure_true)) / 1000);

fprintf('换算为 MPa 后 RMS：%.6e\n', ...
    pressure_rms_raw / 1000);

fprintf('\n假设原始数据单位为 Pa：\n');
fprintf('换算为 MPa 后 RMS：%.6e\n', ...
    pressure_rms_raw / 1e6);

fprintf('\n位移模长拟合出的等效压力：%.6e\n', ...
    p_magnitude);

fprintf('========================================\n');