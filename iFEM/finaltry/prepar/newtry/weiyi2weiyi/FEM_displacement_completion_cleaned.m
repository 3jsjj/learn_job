% =========================================================================
% 物理约束有限元位移补全
%
% 已知：少量测点的 U1/U2/U3、固定边界位移为 0、有限元刚度矩阵 K。
% 未知：其余节点的平移与转动自由度。
%
% 求解：
%     min  1/2 * U * K * U
%     s.t. U(measured translation DOFs) = U_measured
%          U(fixed DOFs) = 0
%
% 分块后求解：K_uu * U_u = -K_uk * U_k
% =========================================================================

clear;
clc;

%% 1. 用户设置

inp_filepath = 'Zhimian_addnodes_cut.inp';
mtx_filepath = 'Job-5-1_STIF2.mtx';
measured_filepath = 'Abaqus_Nodal_U_and_Pressure.csv';
fixed_filepath = 'fixed_nodes.csv';
output_filepath = 'FEM_Physics_Constrained_Displacement_Completion.csv';

% INP/HyperMesh 坐标 -> Abaqus ODB 全局坐标。
R_inp_to_odb = diag([1, -1, -1]);

num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];
fixed_dof_components = 1:6;

% K_uu 的小对角稳定项。
matrix_regularization_factor = 1e-12;

%% 2. 读取 INP 表面网格

[node_ids, node_coords, exterior_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath);

node_coords = node_coords * R_inp_to_odb;

surface_node_rows = unique(exterior_triangles(:));
surface_node_ids = node_ids(surface_node_rows);
surface_coords = node_coords(surface_node_rows, 1:3);
N_surface_nodes = numel(surface_node_rows);

node_row_to_surface_row = zeros(size(node_coords, 1), 1);
node_row_to_surface_row(surface_node_rows) = ...
    (1:N_surface_nodes)';

surface_triangles_local = ...
    node_row_to_surface_row(exterior_triangles);

surface_triangles_local = surface_triangles_local( ...
    all(surface_triangles_local > 0, 2), :);

if isempty(surface_triangles_local)
    error('未找到可用的外表面三角形。');
end

fprintf('读取模型节点：%d\n', numel(node_ids));
fprintf('外表面节点：%d\n', N_surface_nodes);
fprintf('外表面三角形：%d\n', size(surface_triangles_local, 1));

%% 3. 读取并组装 Abaqus 刚度矩阵

try
    mtx_data = readmatrix(mtx_filepath, 'FileType', 'text');
catch
    mtx_data = load(mtx_filepath);
end

mtx_data = mtx_data(all(isfinite(mtx_data), 2), :);

if size(mtx_data, 2) < 5
    error('刚度矩阵文件必须至少包含5列：node_i,dof_i,node_j,dof_j,value。');
end

node_i = round(mtx_data(:,1));
dof_i  = round(mtx_data(:,2));
node_j = round(mtx_data(:,3));
dof_j  = round(mtx_data(:,4));
values = mtx_data(:,5);

if any(~ismember(unique([dof_i; dof_j]), 1:num_dofs_per_node))
    error('MTX 中存在超出每节点6自由度范围的自由度编号。');
end

mtx_node_ids = unique([node_i; node_j]);
N_inp_nodes = numel(node_ids);
node_row_to_mtx_id = nan(N_inp_nodes, 1);

if all(ismember(mtx_node_ids, node_ids))
    % MTX 直接使用 INP NodeID。
    [~, mapped_rows] = ismember(mtx_node_ids, node_ids);
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = INP NodeID';

elseif all(mtx_node_ids >= 1) && ...
       all(mtx_node_ids <= N_inp_nodes) && ...
       all(mtx_node_ids == round(mtx_node_ids))
    % MTX 使用 node_ids/node_coords 的行号。
    mapped_rows = mtx_node_ids;
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = INP节点行号';

else
    error(['无法识别 MTX 节点编号与 INP 节点之间的关系。', newline, ...
           'MTX NodeID 既不是 INP NodeID，也不是有效节点行号。']);
end

row_idx = num_dofs_per_node .* (node_i - 1) + dof_i;
col_idx = num_dofs_per_node .* (node_j - 1) + dof_j;
max_dof = max([row_idx; col_idx]);

K_input = sparse(row_idx, col_idx, values, max_dof, max_dof);

has_upper = nnz(triu(K_input, 1)) > 0;
has_lower = nnz(tril(K_input, -1)) > 0;

if xor(has_upper, has_lower)
    K_raw = K_input + K_input' ...
        - spdiags(diag(K_input), 0, max_dof, max_dof);
else
    symmetry_error = norm(K_input - K_input', 'fro') / ...
        max(norm(K_input, 'fro'), eps);

    if symmetry_error > 1e-8
        warning('MTX 不是严格对称矩阵，当前使用 (K+K'')/2 对称化。');
    end

    K_raw = (K_input + K_input') / 2;
end

[nonzero_row, nonzero_col] = find(K_raw);
existing_dofs = unique([nonzero_row; nonzero_col]);

fprintf('\n========== 刚度矩阵信息 ==========\n');
fprintf('节点映射模式：%s\n', node_mapping_mode);
fprintf('MTX 节点数量：%d\n', numel(mtx_node_ids));
fprintf('原始矩阵尺寸：%d x %d\n', size(K_raw,1), size(K_raw,2));
fprintf('实际存在的自由度数量：%d\n', numel(existing_dofs));
fprintf('矩阵非零项数量：%d\n', nnz(K_raw));
fprintf('==================================\n');

%% 4. 读取固定节点并建立活动自由度

fixed_data = readmatrix(fixed_filepath);
fixed_nodes_input = unique(round(fixed_data(isfinite(fixed_data))));
fixed_nodes_input = fixed_nodes_input(fixed_nodes_input > 0);
fixed_nodes_input = fixed_nodes_input(:);

[fixed_found_in_inp, fixed_node_rows] = ...
    ismember(fixed_nodes_input, node_ids);

if any(~fixed_found_in_inp)
    error('%d 个 fixed_nodes.csv 节点不在 INP 中。', ...
        nnz(~fixed_found_in_inp));
end

fixed_nodes_mtx = node_row_to_mtx_id(fixed_node_rows);

if any(~isfinite(fixed_nodes_mtx))
    error('部分固定节点没有 MTX 节点映射。');
end

fixed_nodes_mtx = round(fixed_nodes_mtx);

fixed_global_dofs_matrix = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* (fixed_nodes_mtx(:) - 1), ...
    fixed_dof_components);

fixed_global_dofs = fixed_global_dofs_matrix(:);
constrained_dofs_existing = intersect( ...
    fixed_global_dofs, existing_dofs);

active_dofs = setdiff(existing_dofs, constrained_dofs_existing);
K_active = K_raw(active_dofs, active_dofs);
N_active_dof = numel(active_dofs);

fprintf('\n========== 固定边界信息 ==========\n');
fprintf('固定节点数量：%d\n', numel(fixed_nodes_input));
fprintf('理论固定自由度：%d\n', ...
    numel(fixed_nodes_input) * numel(fixed_dof_components));
fprintf('实际删除固定自由度：%d\n', ...
    numel(constrained_dofs_existing));
fprintf('剩余活动自由度：%d\n', N_active_dof);
fprintf('==================================\n');

%% 5. 读取测点 U1/U2/U3

measured_data = readmatrix(measured_filepath);

if size(measured_data, 2) < 4
    error('测点文件至少需要四列：NodeID,U1,U2,U3。');
end

measured_ids_raw = measured_data(:,1);
U_measured_raw = measured_data(:,2:4);

valid_measured_rows = ...
    isfinite(measured_ids_raw) & ...
    all(isfinite(U_measured_raw), 2);

measured_ids_raw = round(measured_ids_raw(valid_measured_rows));
U_measured_raw = U_measured_raw(valid_measured_rows, :);

% 同一 NodeID 重复时对位移取平均值。
[measured_ids_unique, ~, measured_group] = unique(measured_ids_raw);
U_measured_unique = zeros(numel(measured_ids_unique), 3);

for component_index = 1:3
    U_measured_unique(:,component_index) = accumarray( ...
        measured_group, ...
        U_measured_raw(:,component_index), ...
        [numel(measured_ids_unique),1], ...
        @mean);
end

[measured_found_in_inp, measured_node_rows] = ...
    ismember(measured_ids_unique, node_ids);

if any(~measured_found_in_inp)
    warning('%d 个测点 NodeID 不在 INP 中，已剔除。', ...
        nnz(~measured_found_in_inp));
end

measured_ids = measured_ids_unique(measured_found_in_inp);
measured_node_rows = measured_node_rows(measured_found_in_inp);
U_measured_mat = U_measured_unique(measured_found_in_inp, :);

measured_mtx_ids = node_row_to_mtx_id(measured_node_rows);
measured_has_mtx_mapping = isfinite(measured_mtx_ids);

measured_global_dofs = nan(numel(measured_ids), 3);
measured_global_dofs(measured_has_mtx_mapping,:) = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* ...
        (measured_mtx_ids(measured_has_mtx_mapping) - 1), ...
    translation_dofs);

% 三个平移自由度必须全部属于活动自由度。
valid_measured_active = ...
    measured_has_mtx_mapping & ...
    all(ismember(measured_global_dofs, active_dofs), 2);

if any(~valid_measured_active)
    warning(['%d 个测点没有完整活动平移自由度，', ...
             '可能属于固定边界或缺少 MTX 映射，已剔除。'], ...
        nnz(~valid_measured_active));
end

measured_ids = measured_ids(valid_measured_active);
measured_node_rows = measured_node_rows(valid_measured_active);
measured_mtx_ids = measured_mtx_ids(valid_measured_active);
measured_global_dofs = measured_global_dofs(valid_measured_active,:);
U_measured_mat = U_measured_mat(valid_measured_active,:);

if isempty(measured_ids)
    error('过滤后没有有效测点。');
end

[~, measured_local_dofs_matrix] = ...
    ismember(measured_global_dofs, active_dofs);

measured_local_dofs = reshape( ...
    measured_local_dofs_matrix', [], 1);

U_measured_vector = reshape(U_measured_mat', [], 1);

% 防止同一活动自由度重复约束。
[known_local_dofs, unique_known_index] = ...
    unique(measured_local_dofs, 'stable');
U_known = U_measured_vector(unique_known_index);

if numel(known_local_dofs) ~= numel(measured_local_dofs)
    warning('测点中存在重复自由度约束，当前保留第一次出现的值。');
end

fprintf('\n========== 测点位移信息 ==========\n');
fprintf('有效测点节点数量：%d\n', numel(measured_ids));
fprintf('严格约束平移自由度：%d\n', numel(known_local_dofs));
fprintf('输入测点位移范数：%.6e\n', norm(U_known));
fprintf('==================================\n');

%% 6. 刚度矩阵最小应变能位移补全

all_active_local_dofs = (1:N_active_dof)';
unknown_local_dofs = setdiff( ...
    all_active_local_dofs, known_local_dofs);

U_active_completed = zeros(N_active_dof, 1);
U_active_completed(known_local_dofs) = U_known;

if isempty(unknown_local_dofs)
    warning('所有活动自由度均已知，无需补全。');
else
    K_uu = K_active(unknown_local_dofs, unknown_local_dofs);
    K_uk = K_active(unknown_local_dofs, known_local_dofs);

    rhs_unknown = -K_uk * U_known;

    diagonal_values = abs(full(diag(K_uu)));
    positive_diagonal = diagonal_values(diagonal_values > 0 & ...
        isfinite(diagonal_values));

    if isempty(positive_diagonal)
        diagonal_scale = 1;
    else
        diagonal_scale = median(positive_diagonal);
    end

    regularization_value = ...
        matrix_regularization_factor * diagonal_scale;

    K_uu_regularized = K_uu + ...
        regularization_value * speye(size(K_uu,1));

    fprintf('\n正在求解未知自由度：%d ...\n', ...
        numel(unknown_local_dofs));
    fprintf('对角稳定项：%.6e\n', regularization_value);

    U_unknown = K_uu_regularized \ rhs_unknown;
    U_active_completed(unknown_local_dofs) = U_unknown;
end

% 再次强制写回测量值，消除数值舍入误差。
U_active_completed(known_local_dofs) = U_known;

fprintf('位移补全完成：已知自由度 %d，补全自由度 %d。\n', ...
    numel(known_local_dofs), numel(unknown_local_dofs));

%% 7. 提取外表面节点位移

surface_mtx_ids = node_row_to_mtx_id(surface_node_rows);
U_surface_6D = nan(N_surface_nodes, 6);

for surface_index = 1:N_surface_nodes

    current_mtx_id = surface_mtx_ids(surface_index);

    if ~isfinite(current_mtx_id)
        continue;
    end

    current_global_dofs = ...
        num_dofs_per_node * (current_mtx_id - 1) + (1:6);

    for component_index = 1:6
        current_global_dof = current_global_dofs(component_index);

        [is_active, active_location] = ...
            ismember(current_global_dof, active_dofs);

        if is_active
            U_surface_6D(surface_index, component_index) = ...
                U_active_completed(active_location);

        elseif ismember(current_global_dof, ...
                constrained_dofs_existing)
            U_surface_6D(surface_index, component_index) = 0;
        end
    end
end

U_surface_3D = U_surface_6D(:,1:3);
U_surface_rotation = U_surface_6D(:,4:6);
U_magnitude = vecnorm(U_surface_3D, 2, 2);

%% 8. 输出 CSV

result_table = table( ...
    surface_node_ids(:), ...
    surface_coords(:,1), ...
    surface_coords(:,2), ...
    surface_coords(:,3), ...
    U_surface_3D(:,1), ...
    U_surface_3D(:,2), ...
    U_surface_3D(:,3), ...
    U_surface_rotation(:,1), ...
    U_surface_rotation(:,2), ...
    U_surface_rotation(:,3), ...
    U_magnitude, ...
    'VariableNames', { ...
        'NodeID', 'X', 'Y', 'Z', ...
        'U1_Completed', 'U2_Completed', 'U3_Completed', ...
        'UR1_Completed', 'UR2_Completed', 'UR3_Completed', ...
        'U_Magnitude'});

writetable(result_table, output_filepath);
fprintf('结果已保存：%s\n', output_filepath);
fprintf('程序运行完成。\n');
