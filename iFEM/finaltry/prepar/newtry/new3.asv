% =========================================================================
% 物理约束有限元位移补全（不反演压力、不输入载荷）
%
% 核心思想：
%   已知：少量测点的 U1/U2/U3、固定边界位移为 0、有限元刚度矩阵 K。
%   未知：其余节点的平移与转动自由度。
%
%   在严格满足已知位移的条件下，求使结构总应变能最小的位移场：
%
%       min  1/2 * U' * K * U
%       s.t. U(measured translation DOFs) = U_measured
%            U(fixed DOFs) = 0
%
%   将活动自由度分成已知 k 和未知 u 后：
%
%       K_uu * U_u = -K_uk * U_k
%
% 注意：
%   1) 本程序不读取压力、不构造法向/面积、不构造载荷矩阵 G；
%   2) 测点位移被作为严格位移约束，因此会原样保留；
%   3) 该解等价于“未知自由度没有直接外力时的最低应变能解”；
%      若真实载荷分布在大量未知节点上，结果仍可能存在误差；
%   4) 当前壳模型按每节点 6 自由度处理：U1,U2,U3,UR1,UR2,UR3。
% =========================================================================

clear;
clc;

%% 1. 用户设置

inp_filepath = 'Zhimian_addnodes_cut.inp';
mtx_filepath = 'Job-5-1_STIF2.mtx';

% 至少四列：NodeID,U1,U2,U3。
% 后面的压力或其他列会被忽略。
measured_filepath = 'Abaqus_Nodal_U_and_Pressure.csv';

% 固定节点 NodeID；默认约束每个固定节点的 1:6 自由度。
fixed_filepath = 'fixed_nodes.csv';

% 可选：完整 Abaqus 表面位移，用于验证补全结果。
% 至少四列：NodeID,U1,U2,U3。
reference_filepath = 'Abaqus_All_Surface_Displacement.csv';
compare_with_reference_if_available = true;

% INP/HyperMesh 坐标 -> Abaqus ODB 全局坐标。
% 仅用于坐标和绘图；若位移直接来自 ODB，不要再变换 U1/U2/U3。
R_inp_to_odb = diag([1, -1, -1]);

num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];
rotation_dofs = [4, 5, 6];
fixed_dof_components = 1:6;

% 对 K_uu 添加极小对角稳定项：
%   K_uu_regularized = K_uu + factor * diagonal_scale * I
% 正常情况下保持很小；若求解提示奇异，可逐步增大到 1e-10 或 1e-8。
matrix_regularization_factor = 1e-12;

% 位移显示单位。
displacement_unit = 'mm';

save_figures = true;

%% 2. 读取 INP 表面网格

% 该函数沿用你原程序中的网格读取函数。
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

% K*U 可解释为维持该补全位移场所需要的等效节点力。
% 本程序不把它当作真实外载荷，只用于诊断。
F_inferred_active = K_active * U_active_completed;

if isempty(unknown_local_dofs)
    unknown_equilibrium_residual = 0;
else
    unknown_force_residual = F_inferred_active(unknown_local_dofs);
    unknown_equilibrium_residual = ...
        norm(unknown_force_residual) / ...
        max(norm(F_inferred_active), eps);
end

measured_reproduction_error = ...
    norm(U_active_completed(known_local_dofs) - U_known) / ...
    max(norm(U_known), eps);

strain_energy = 0.5 * real( ...
    U_active_completed' * K_active * U_active_completed);

fprintf('\n========== 位移补全诊断 ==========\n');
fprintf('已知活动自由度：%d\n', numel(known_local_dofs));
fprintf('补全未知自由度：%d\n', numel(unknown_local_dofs));
fprintf('测点位移重现相对误差：%.6e\n', ...
    measured_reproduction_error);
fprintf('未知自由度平衡相对残差：%.6e\n', ...
    unknown_equilibrium_residual);
fprintf('补全位移场应变能：%.6e\n', strain_energy);
fprintf('==================================\n');

%% 7. 提取全部外表面节点的 U1/U2/U3/UR1/UR2/UR3

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

% 标记测点和固定节点。
[is_measured_surface, ~] = ...
    ismember(surface_node_ids, measured_ids);
[is_fixed_surface, ~] = ...
    ismember(surface_node_ids, fixed_nodes_input);

% 检查测点值在表面输出中是否被精确保留。
[measured_on_surface, measured_surface_location] = ...
    ismember(measured_ids, surface_node_ids);

if any(measured_on_surface)
    measured_output_error = ...
        U_surface_3D(measured_surface_location(measured_on_surface),:) - ...
        U_measured_mat(measured_on_surface,:);

    measured_surface_relative_error = ...
        norm(measured_output_error, 'fro') / ...
        max(norm(U_measured_mat(measured_on_surface,:), 'fro'), eps);
else
    measured_surface_relative_error = NaN;
end

fprintf('表面测点位移保持误差：%.6e\n', ...
    measured_surface_relative_error);

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
    is_measured_surface(:), ...
    is_fixed_surface(:), ...
    'VariableNames', { ...
        'NodeID', ...
        'X','Y','Z', ...
        'U1_Completed','U2_Completed','U3_Completed', ...
        'UR1_Completed','UR2_Completed','UR3_Completed', ...
        'U_Magnitude', ...
        'IsMeasuredNode','IsFixedNode'});

output_filepath = ...
    'FEM_Physics_Constrained_Displacement_Completion.csv';
writetable(result_table, output_filepath);

fprintf('结果已保存：%s\n', output_filepath);

%% 9. 绘制完整外表面位移云图

plot_values = { ...
    U_surface_3D(:,1), ...
    U_surface_3D(:,2), ...
    U_surface_3D(:,3), ...
    U_magnitude};

plot_names = {'U1','U2','U3','UMagnitude'};
plot_titles = { ...
    '物理约束有限元位移补全 U1', ...
    '物理约束有限元位移补全 U2', ...
    '物理约束有限元位移补全 U3', ...
    '物理约束有限元位移补全模长'};

for plot_index = 1:numel(plot_values)

    current_value = plot_values{plot_index};

    figure( ...
        'Name', plot_titles{plot_index}, ...
        'Color', 'w', ...
        'Position', [100,100,900,700]);

    trisurf( ...
        surface_triangles_local, ...
        surface_coords(:,1), ...
        surface_coords(:,2), ...
        surface_coords(:,3), ...
        current_value, ...
        'FaceColor','interp', ...
        'EdgeColor','none');

    hold on;

    scatter3( ...
        surface_coords(is_measured_surface,1), ...
        surface_coords(is_measured_surface,2), ...
        surface_coords(is_measured_surface,3), ...
        34, 'k', 'filled');

    scatter3( ...
        surface_coords(is_fixed_surface,1), ...
        surface_coords(is_fixed_surface,2), ...
        surface_coords(is_fixed_surface,3), ...
        20, 'r', 'filled');

    axis equal;
    axis tight;
    grid on;
    box on;
    view(45,30);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(plot_titles{plot_index});

    cb = colorbar;
    ylabel(cb, sprintf('%s (%s)', ...
        plot_names{plot_index}, displacement_unit));

    colormap(jet(256));

    if plot_index <= 3
        finite_values = current_value(isfinite(current_value));
        if ~isempty(finite_values)
            current_limit = max(abs(finite_values));
            if current_limit > 0
                caxis([-current_limit,current_limit]);
            end
        end
    end

    legend('补全曲面','测量节点','固定节点', ...
        'Location','best');
    hold off;

    if save_figures
        output_image_name = sprintf( ...
            'FEM_Completed_%s_Cloud.png', ...
            plot_names{plot_index});

        exportgraphics(gcf, output_image_name, ...
            'Resolution',300);
        fprintf('已保存：%s\n', output_image_name);
    end
end

%% 10. 可选：与完整 Abaqus 表面位移比较

if compare_with_reference_if_available && ...
        isfile(reference_filepath)

    reference_data = readmatrix(reference_filepath);

    if size(reference_data,2) < 4
        warning('完整位移参考文件少于四列，跳过比较。');
    else
        reference_ids_raw = reference_data(:,1);
        U_reference_raw = reference_data(:,2:4);

        valid_reference_rows = ...
            isfinite(reference_ids_raw) & ...
            all(isfinite(U_reference_raw),2);

        reference_ids = round( ...
            reference_ids_raw(valid_reference_rows));
        U_reference = U_reference_raw(valid_reference_rows,:);

        [matched_surface, reference_location] = ...
            ismember(surface_node_ids, reference_ids);

        valid_completed_surface = ...
            all(isfinite(U_surface_3D),2);

        compare_mask = ...
            matched_surface & valid_completed_surface;
        compare_rows = find(compare_mask);

        U_completed_compare = ...
            U_surface_3D(compare_rows,:);
        U_reference_compare = U_reference( ...
            reference_location(compare_rows),:);

        completion_error = ...
            U_completed_compare - U_reference_compare;

        fprintf('\n');
        fprintf('========== 完整位移场比较 ==========\n');
        fprintf('匹配节点数量：%d / %d\n', ...
            numel(compare_rows), N_surface_nodes);

        component_names = {'U1','U2','U3'};

        for component_index = 1:3
            relative_error = ...
                norm(completion_error(:,component_index)) / ...
                max(norm(U_reference_compare(:,component_index)),eps);

            if std(U_completed_compare(:,component_index)) > eps && ...
                    std(U_reference_compare(:,component_index)) > eps
                current_correlation = corrcoef( ...
                    U_completed_compare(:,component_index), ...
                    U_reference_compare(:,component_index));
                correlation_value = current_correlation(1,2);
            else
                correlation_value = NaN;
            end

            fprintf('%s 相对L2误差：%.2f %%\n', ...
                component_names{component_index}, ...
                100 * relative_error);
            fprintf('%s 相关系数：%.6f\n', ...
                component_names{component_index}, ...
                correlation_value);
        end

        overall_relative_error = ...
            norm(completion_error,'fro') / ...
            max(norm(U_reference_compare,'fro'),eps);

        overall_cosine_similarity = ...
            dot(U_completed_compare(:), ...
                U_reference_compare(:)) / ...
            max(norm(U_completed_compare(:)) * ...
                norm(U_reference_compare(:)),eps);

        best_scale = ...
            dot(U_completed_compare(:), ...
                U_reference_compare(:)) / ...
            max(dot(U_completed_compare(:), ...
                    U_completed_compare(:)),eps);

        scaled_completed = best_scale * U_completed_compare;
        scaled_relative_error = ...
            norm(scaled_completed - U_reference_compare,'fro') / ...
            max(norm(U_reference_compare,'fro'),eps);

        fprintf('整体相对L2误差：%.2f %%\n', ...
            100 * overall_relative_error);
        fprintf('整体余弦相似度：%.6f\n', ...
            overall_cosine_similarity);
        fprintf('补全场到Abaqus最佳比例：%.6e\n', ...
            best_scale);
        fprintf('消除整体比例后的误差：%.2f %%\n', ...
            100 * scaled_relative_error);
        fprintf('补全位移场范数：%.6e\n', ...
            norm(U_completed_compare,'fro'));
        fprintf('Abaqus参考场范数：%.6e\n', ...
            norm(U_reference_compare,'fro'));
        fprintf('====================================\n');

        error_magnitude = vecnorm(completion_error,2,2);

        comparison_table = table( ...
            surface_node_ids(compare_rows), ...
            U_reference_compare(:,1), ...
            U_completed_compare(:,1), ...
            completion_error(:,1), ...
            U_reference_compare(:,2), ...
            U_completed_compare(:,2), ...
            completion_error(:,2), ...
            U_reference_compare(:,3), ...
            U_completed_compare(:,3), ...
            completion_error(:,3), ...
            error_magnitude, ...
            is_measured_surface(compare_rows), ...
            is_fixed_surface(compare_rows), ...
            'VariableNames', { ...
                'NodeID', ...
                'U1_Abaqus','U1_Completed','U1_Error', ...
                'U2_Abaqus','U2_Completed','U2_Error', ...
                'U3_Abaqus','U3_Completed','U3_Error', ...
                'ErrorMagnitude', ...
                'IsMeasuredNode','IsFixedNode'});

        writetable(comparison_table, ...
            'Abaqus_vs_FEM_Completed_Displacement.csv');

        [~, worst_order] = sort(error_magnitude,'descend');
        worst_count = min(10,numel(worst_order));

        if worst_count > 0
            fprintf('\n误差最大的节点：\n');
            disp(comparison_table( ...
                worst_order(1:worst_count), ...
                {'NodeID','ErrorMagnitude', ...
                 'IsMeasuredNode','IsFixedNode'}));
        end

        fprintf('比较结果已保存：');
        fprintf('Abaqus_vs_FEM_Completed_Displacement.csv\n');

        % 误差云图。
        error_surface = nan(N_surface_nodes,1);
        error_surface(compare_rows) = error_magnitude;

        figure( ...
            'Name','FEM位移补全误差模长', ...
            'Color','w', ...
            'Position',[100,100,900,700]);

        trisurf( ...
            surface_triangles_local, ...
            surface_coords(:,1), ...
            surface_coords(:,2), ...
            surface_coords(:,3), ...
            error_surface, ...
            'FaceColor','interp', ...
            'EdgeColor','none');

        axis equal;
        axis tight;
        grid on;
        box on;
        view(45,30);
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('物理约束有限元位移补全误差模长');
        cb = colorbar;
        ylabel(cb, ['||U_{completed}-U_{Abaqus}|| (', ...
            displacement_unit, ')']);
        colormap(jet(256));

        if save_figures
            exportgraphics(gcf, ...
                'FEM_Completion_Error_Cloud.png', ...
                'Resolution',300);
        end
    end
end

fprintf('\n程序运行完成。\n');
