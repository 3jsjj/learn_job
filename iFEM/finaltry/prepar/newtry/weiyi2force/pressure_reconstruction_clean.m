% =========================================================================
% 完整六自由度位移 -> 等效节点力 -> 自动识别外表面压力（精简主线版）
%
% 主流程：
%   1) 读取 INP 外表面网格、MTX 刚度矩阵、固定节点和完整六自由度位移；
%   2) 计算等效节点广义力 F = K * U；
%   3) 在全部外表面建立一致压力载荷矩阵 G；
%   4) 通过非负、平滑、稀疏正则化反演节点压力；
%   5) 自动比较正/反法向并输出最终压力结果。
%
% 输出：
%   Pressure_Reconstruction.csv
%   Pressure_Reconstruction.png（save_figure = true 时）
%
% 依赖：
%   read_abaqus_surface_mesh.m
%   返回 [node_ids, node_coords, exterior_triangles]
%   exterior_triangles 中保存节点行号，而不是 NodeID。
% =========================================================================

clear;
clc;
close all;

%% 1. 用户设置

inp_filepath = 'Zhimian_addnodes_cut.inp';
mtx_filepath = 'Job-5-1_STIF2.mtx';
displacement_filepath = 'Abaqus_All_Nodal_U_UR.csv';
fixed_filepath = 'fixed_nodes.csv';

% INP/HyperMesh 坐标 -> Abaqus ODB 全局坐标。
% ODB 导出的位移和转角不再做此变换。
R_inp_to_odb = diag([1, -1, -1]);

num_dofs_per_node = 6;
translation_dofs = 1:3;
% 压力方向：-1 为沿外法向反方向，+1 为沿外法向。
auto_choose_pressure_sign = true;
manual_pressure_sign = -1;
pressure_sign_candidates = [-1, +1];
orient_surface_normals_outward = true;

% 正则化参数（G 和 D 会自动归一化）。
smoothness_weight = 2e-2;
sparsity_fraction = 5e-3;
ridge_weight = 1e-10;

% FISTA 设置。
solver_max_iterations = 6000;
solver_tolerance = 1e-9;

% 仅用于标记活动受压节点，不影响求解。
active_pressure_relative_threshold = 0.02;

% false：位移缺失时报错；true：缺失自由度按 0 处理。
allow_missing_displacement_as_zero = false;

pressure_unit = 'MPa';
save_figure = true;

%% 2. 读取外表面网格

[node_ids, node_coords, exterior_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath);

node_ids = node_ids(:);
node_coords = node_coords(:, 1:3) * R_inp_to_odb;
N_model_nodes = numel(node_ids);

if isempty(exterior_triangles)
    error('没有读取到外表面三角形。');
end

%% 3. 读取并组装刚度矩阵

try
    mtx_data = readmatrix(mtx_filepath, 'FileType', 'text');
catch
    mtx_data = load(mtx_filepath);
end

mtx_data = mtx_data(all(isfinite(mtx_data), 2), :);

if size(mtx_data, 2) < 5
    error('MTX 文件至少应包含 5 列：node_i,dof_i,node_j,dof_j,value。');
end

node_i = round(mtx_data(:, 1));
dof_i = round(mtx_data(:, 2));
node_j = round(mtx_data(:, 3));
dof_j = round(mtx_data(:, 4));
values = mtx_data(:, 5);

matrix_dof_components = unique([dof_i; dof_j]);
if any(~ismember(matrix_dof_components, 1:num_dofs_per_node))
    error('MTX 中出现了超出每节点 6 自由度范围的编号。');
end

mtx_node_ids = unique([node_i; node_j]);
node_row_to_mtx_id = nan(N_model_nodes, 1);

if all(ismember(mtx_node_ids, node_ids))
    [~, mapped_rows] = ismember(mtx_node_ids, node_ids);
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
elseif all(mtx_node_ids >= 1) && ...
       all(mtx_node_ids <= N_model_nodes) && ...
       all(mtx_node_ids == round(mtx_node_ids))
    mapped_rows = mtx_node_ids;
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
else
    error(['无法识别 MTX NodeID 与 INP 节点的关系。', newline, ...
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
        warning('输入刚度矩阵不是严格对称矩阵，已进行对称化。');
    end
    K_raw = (K_input + K_input') / 2;
end

[nonzero_row, nonzero_col] = find(K_raw);
existing_dofs = unique([nonzero_row; nonzero_col]);

%% 4. 读取固定节点

fixed_data = readmatrix(fixed_filepath);
fixed_node_ids = unique(round(fixed_data(isfinite(fixed_data))));
fixed_node_ids = fixed_node_ids(fixed_node_ids > 0);
fixed_node_ids = fixed_node_ids(:);

[fixed_found, fixed_node_rows] = ismember(fixed_node_ids, node_ids);
if any(~fixed_found)
    error('%d 个固定节点不在 INP 节点中。', nnz(~fixed_found));
end

fixed_nodes_mtx = node_row_to_mtx_id(fixed_node_rows);
if any(~isfinite(fixed_nodes_mtx))
    error('部分固定节点没有 MTX 映射。');
end

%% 5. 读取完整六自由度位移并组装 U

displacement_data = readmatrix(displacement_filepath);

if size(displacement_data, 2) < 7
    error('位移文件必须至少包含：NodeID,U1,U2,U3,UR1,UR2,UR3。');
end

displacement_ids_raw = displacement_data(:, 1);
U6_raw = displacement_data(:, 2:7);
valid_id_rows = isfinite(displacement_ids_raw);

displacement_ids_raw = round(displacement_ids_raw(valid_id_rows));
U6_raw = U6_raw(valid_id_rows, :);

% 合并重复 NodeID：每个有限分量取平均值。
[displacement_ids_unique, ~, displacement_group] = ...
    unique(displacement_ids_raw);
U6_unique = nan(numel(displacement_ids_unique), 6);

for component_index = 1:6
    for group_index = 1:numel(displacement_ids_unique)
        current_values = U6_raw( ...
            displacement_group == group_index, component_index);
        current_values = current_values(isfinite(current_values));
        if ~isempty(current_values)
            U6_unique(group_index, component_index) = mean(current_values);
        end
    end
end

[displacement_found, displacement_node_rows] = ...
    ismember(displacement_ids_unique, node_ids);

if any(~displacement_found)
    warning('%d 个位移 NodeID 不在 INP 中，已忽略。', ...
        nnz(~displacement_found));
end

displacement_node_rows = displacement_node_rows(displacement_found);
U6_data = U6_unique(displacement_found, :);

U_full = zeros(max_dof, 1);
U_dof_has_value = false(max_dof, 1);

for data_index = 1:numel(displacement_node_rows)
    node_row = displacement_node_rows(data_index);
    mtx_node_id = node_row_to_mtx_id(node_row);

    if ~isfinite(mtx_node_id)
        continue;
    end

    absolute_dofs = ...
        num_dofs_per_node * (mtx_node_id - 1) + (1:6);

    for component_index = 1:6
        current_dof = absolute_dofs(component_index);
        if current_dof > max_dof
            continue;
        end

        current_value = U6_data(data_index, component_index);
        if isfinite(current_value)
            U_full(current_dof) = current_value;
            U_dof_has_value(current_dof) = true;
        end
    end
end

missing_existing_dofs = existing_dofs(~U_dof_has_value(existing_dofs));
if ~isempty(missing_existing_dofs)
    if allow_missing_displacement_as_zero
        warning('%d 个 MTX 自由度缺少位移，已按 0 处理。', ...
            numel(missing_existing_dofs));
    else
        error(['完整位移文件没有覆盖全部 MTX 自由度。', newline, ...
               '请确认 ODB 已导出 U1,U2,U3,UR1,UR2,UR3。']);
    end
end

%% 6. 计算等效节点广义力 F = K * U

F_equivalent_full = K_raw * U_full;

%% 7. 建立全部外表面的压力未知量

pressure_triangles = exterior_triangles;
if orient_surface_normals_outward
    pressure_triangles = orient_triangles_by_model_centroid( ...
        pressure_triangles, node_coords);
end

pressure_node_rows = unique(pressure_triangles(:));
pressure_node_ids = node_ids(pressure_node_rows);
pressure_coords = node_coords(pressure_node_rows, :);
N_pressure_nodes = numel(pressure_node_rows);
N_pressure_triangles = size(pressure_triangles, 1);

pressure_node_mtx_ids = node_row_to_mtx_id(pressure_node_rows);
if any(~isfinite(pressure_node_mtx_ids))
    error('部分外表面节点没有 MTX 映射。');
end

pressure_node_is_fixed = ismember(pressure_node_rows, fixed_node_rows);

%% 8. 建立一致压力载荷矩阵 G

pressure_node_index_of_model_row = zeros(N_model_nodes, 1);
pressure_node_index_of_model_row(pressure_node_rows) = ...
    (1:N_pressure_nodes)';

estimated_entries = 27 * N_pressure_triangles;
G_rows = zeros(estimated_entries, 1);
G_cols = zeros(estimated_entries, 1);
G_vals = zeros(estimated_entries, 1);
entry_cursor = 0;

consistent_triangle_weights = ...
    [2, 1, 1; 1, 2, 1; 1, 1, 2] / 12;

for triangle_index = 1:N_pressure_triangles
    tri_rows = pressure_triangles(triangle_index, :);
    tri_coords = node_coords(tri_rows, :);

    area_vector = 0.5 * cross( ...
        tri_coords(2, :) - tri_coords(1, :), ...
        tri_coords(3, :) - tri_coords(1, :));

    if norm(area_vector) <= eps
        error('外表面第 %d 个三角形面积为 0。', triangle_index);
    end

    local_pressure_indices = ...
        pressure_node_index_of_model_row(tri_rows);

    for local_force_node = 1:3
        force_node_row = tri_rows(local_force_node);
        force_mtx_id = node_row_to_mtx_id(force_node_row);
        force_absolute_dofs = ...
            num_dofs_per_node * (force_mtx_id - 1) + translation_dofs;

        for local_pressure_node = 1:3
            pressure_column = ...
                local_pressure_indices(local_pressure_node);
            coefficient_vector = ...
                consistent_triangle_weights( ...
                    local_force_node, local_pressure_node) * area_vector;

            current_entries = entry_cursor + (1:3);
            G_rows(current_entries) = force_absolute_dofs(:);
            G_cols(current_entries) = pressure_column;
            G_vals(current_entries) = coefficient_vector(:);
            entry_cursor = entry_cursor + 3;
        end
    end
end

G_outward_full = sparse( ...
    G_rows(1:entry_cursor), ...
    G_cols(1:entry_cursor), ...
    G_vals(1:entry_cursor), ...
    max_dof, N_pressure_nodes);

[pressure_node_normals, pressure_node_areas] = ...
    compute_pressure_node_geometry( ...
        pressure_node_rows, pressure_triangles, node_coords);

%% 9. 选择用于压力拟合的自由表面平移自由度

pressure_translation_dofs_matrix = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* (pressure_node_mtx_ids(:) - 1), ...
    translation_dofs);

fit_pressure_node_mask = ...
    ~pressure_node_is_fixed & ...
    all(ismember(pressure_translation_dofs_matrix, existing_dofs), 2);

force_fit_dofs_matrix = ...
    pressure_translation_dofs_matrix(fit_pressure_node_mask, :);
force_fit_dofs = reshape(force_fit_dofs_matrix', [], 1);

G_outward_fit = G_outward_full(force_fit_dofs, :);
F_fit = F_equivalent_full(force_fit_dofs);

row_has_pressure_sensitivity = ...
    full(sum(abs(G_outward_fit), 2)) > 0;
G_outward_fit = G_outward_fit(row_has_pressure_sensitivity, :);
F_fit = F_fit(row_has_pressure_sensitivity);
force_fit_dofs = force_fit_dofs(row_has_pressure_sensitivity);

if isempty(force_fit_dofs)
    error('没有可用于压力反演的外表面自由节点平移力方程。');
end

%% 10. 构造光滑矩阵并归一化

D_pressure = build_edge_difference_matrix( ...
    pressure_triangles, ...
    pressure_node_index_of_model_row, ...
    N_pressure_nodes);

G_scale = normest(G_outward_fit);
if ~isfinite(G_scale) || G_scale <= eps
    error('压力映射矩阵 G 几乎为零。');
end

G_outward_normalized = G_outward_fit / G_scale;
F_normalized = F_fit / G_scale;

D_scale = normest(D_pressure);
if ~isfinite(D_scale) || D_scale <= eps
    D_normalized = D_pressure;
else
    D_normalized = D_pressure / D_scale;
end

%% 11. 非负稀疏压力反演并选择压力方向

if auto_choose_pressure_sign
    tested_pressure_signs = pressure_sign_candidates(:);
else
    tested_pressure_signs = manual_pressure_sign;
end

N_signs = numel(tested_pressure_signs);
pressure_candidates = cell(N_signs, 1);
solver_info_candidates = cell(N_signs, 1);
sign_relative_residuals = inf(N_signs, 1);
sign_l1_weights = zeros(N_signs, 1);

for sign_index = 1:N_signs
    current_sign = tested_pressure_signs(sign_index);
    G_candidate = current_sign * G_outward_normalized;

    beta_max = max(full(G_candidate' * F_normalized));
    beta_max = max(beta_max, 0);
    current_l1_weight = sparsity_fraction * beta_max;

    [current_pressure, current_info] = ...
        solve_sparse_nonnegative_pressure_fista( ...
            G_candidate, ...
            F_normalized, ...
            D_normalized, ...
            smoothness_weight, ...
            current_l1_weight, ...
            ridge_weight, ...
            solver_max_iterations, ...
            solver_tolerance);

    pressure_candidates{sign_index} = current_pressure;
    solver_info_candidates{sign_index} = current_info;
    sign_relative_residuals(sign_index) = current_info.relative_residual;
    sign_l1_weights(sign_index) = current_l1_weight;
end

[~, best_sign_index] = min(sign_relative_residuals);
pressure_sign_selected = tested_pressure_signs(best_sign_index);
pressure_nodal = pressure_candidates{best_sign_index};
pressure_solver_info = solver_info_candidates{best_sign_index};
sparsity_weight_selected = sign_l1_weights(best_sign_index);
G_full = pressure_sign_selected * G_outward_full;

if ~any(pressure_nodal > 0)
    warning(['压力反演结果全部为 0。请检查位移、刚度矩阵、', ...
             '坐标变换和法向方向，或减小 sparsity_fraction。']);
end

%% 12. 计算拟合质量和活动区域

F_pressure_full = G_full * pressure_nodal;
F_pressure_fit = F_pressure_full(force_fit_dofs);
F_equivalent_fit = F_equivalent_full(force_fit_dofs);
force_fit_residual = F_pressure_fit - F_equivalent_fit;

relative_force_fit_error = ...
    norm(force_fit_residual) / max(norm(F_equivalent_fit), eps);
force_fit_cosine = dot(F_pressure_fit, F_equivalent_fit) / ...
    max(norm(F_pressure_fit) * norm(F_equivalent_fit), eps);

if max(pressure_nodal) > 0
    pressure_relative_to_max = pressure_nodal / max(pressure_nodal);
    is_active_pressure_node = ...
        pressure_relative_to_max >= active_pressure_relative_threshold;
else
    pressure_relative_to_max = zeros(size(pressure_nodal));
    is_active_pressure_node = false(size(pressure_nodal));
end

%% 13. 输出最终压力结果

pressure_result_table = table( ...
    pressure_node_ids(:), ...
    pressure_coords(:, 1), ...
    pressure_coords(:, 2), ...
    pressure_coords(:, 3), ...
    pressure_node_normals(:, 1), ...
    pressure_node_normals(:, 2), ...
    pressure_node_normals(:, 3), ...
    pressure_node_areas(:), ...
    pressure_nodal(:), ...
    pressure_relative_to_max(:), ...
    is_active_pressure_node(:), ...
    pressure_node_is_fixed(:), ...
    'VariableNames', { ...
        'NodeID', 'X', 'Y', 'Z', ...
        'NormalX', 'NormalY', 'NormalZ', ...
        'NodalArea', 'Pressure_Reconstructed', ...
        'Pressure_RelativeToMax', ...
        'IsActivePressureNode', 'IsFixedNode'});

writetable(pressure_result_table, 'Pressure_Reconstruction.csv');

fprintf('\n========== 压力反演完成 ==========\n');
fprintf('压力方向：%+d\n', pressure_sign_selected);
fprintf('候选压力节点：%d\n', N_pressure_nodes);
fprintf('活动压力节点：%d\n', nnz(is_active_pressure_node));
fprintf('最大压力：%.6e %s\n', max(pressure_nodal), pressure_unit);
fprintf('相对拟合误差：%.2f %%\n', 100 * relative_force_fit_error);
fprintf('力向量余弦相似度：%.6f\n', force_fit_cosine);
fprintf('FISTA 迭代次数：%d\n', pressure_solver_info.iterations);
fprintf('稀疏权重 beta：%.6e\n', sparsity_weight_selected);
fprintf('已保存：Pressure_Reconstruction.csv\n');
fprintf('==================================\n');

%% 14. 绘制最终压力云图

pressure_row_to_local = zeros(N_model_nodes, 1);
pressure_row_to_local(pressure_node_rows) = (1:N_pressure_nodes)';
pressure_triangles_local = pressure_row_to_local(pressure_triangles);

figure( ...
    'Name', '自动识别外表面压力', ...
    'Color', 'w', ...
    'Position', [100, 100, 900, 700]);

trisurf( ...
    pressure_triangles_local, ...
    pressure_coords(:, 1), ...
    pressure_coords(:, 2), ...
    pressure_coords(:, 3), ...
    pressure_nodal, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none');

axis equal;
axis tight;
grid on;
box on;
view(45, 30);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('由完整六自由度位移反演的外表面压力');
cb = colorbar;
ylabel(cb, ['Pressure (', pressure_unit, ')']);
colormap(jet(256));

if save_figure
    exportgraphics(gcf, 'Pressure_Reconstruction.png', 'Resolution', 300);
    fprintf('已保存：Pressure_Reconstruction.png\n');
end

%% =========================================================================
% 局部函数
% =========================================================================

function [pressure, info] = ...
    solve_sparse_nonnegative_pressure_fista( ...
        G, F, D, lambda_smooth, lambda_l1, lambda_ridge, ...
        max_iterations, tolerance)
% 求解：
%   min_{p>=0} 1/2||Gp-F||^2
%            + 1/2 lambda_smooth^2 ||Dp||^2
%            + lambda_l1 ||p||_1
%            + 1/2 lambda_ridge ||p||^2

    N_pressure = size(G, 2);
    pressure = zeros(N_pressure, 1);
    extrapolated_pressure = pressure;
    momentum = 1;

    G_lipschitz = normest(G)^2;
    if isempty(D) || nnz(D) == 0
        D_lipschitz = 0;
    else
        D_lipschitz = normest(D)^2;
    end

    lipschitz_constant = ...
        G_lipschitz + ...
        lambda_smooth^2 * D_lipschitz + ...
        lambda_ridge;
    lipschitz_constant = max(lipschitz_constant, eps);

    objective_previous = inf;
    converged = false;

    for iteration = 1:max_iterations
        gradient_value = ...
            G' * (G * extrapolated_pressure - F) + ...
            lambda_smooth^2 * (D' * (D * extrapolated_pressure)) + ...
            lambda_ridge * extrapolated_pressure;

        proximal_argument = ...
            extrapolated_pressure - gradient_value / lipschitz_constant;
        pressure_new = max( ...
            proximal_argument - lambda_l1 / lipschitz_constant, 0);

        momentum_new = (1 + sqrt(1 + 4 * momentum^2)) / 2;
        extrapolated_new = pressure_new + ...
            ((momentum - 1) / momentum_new) * ...
            (pressure_new - pressure);

        current_residual = G * pressure_new - F;
        current_objective = ...
            0.5 * norm(current_residual)^2 + ...
            0.5 * lambda_smooth^2 * norm(D * pressure_new)^2 + ...
            lambda_l1 * sum(pressure_new) + ...
            0.5 * lambda_ridge * norm(pressure_new)^2;

        % 单调重启：目标函数升高时取消动量。
        if current_objective > objective_previous
            extrapolated_new = pressure_new;
            momentum_new = 1;
        end

        relative_change = ...
            norm(pressure_new - pressure) / max(norm(pressure_new), 1);

        pressure = pressure_new;
        extrapolated_pressure = extrapolated_new;
        momentum = momentum_new;
        objective_previous = current_objective;

        if relative_change < tolerance
            converged = true;
            break;
        end
    end

    info.iterations = iteration;
    info.converged = converged;
    info.objective = objective_previous;
    info.relative_residual = ...
        norm(G * pressure - F) / max(norm(F), eps);
    info.lipschitz_constant = lipschitz_constant;
end

function triangles_out = ...
    orient_triangles_by_model_centroid(triangles_in, node_coords)

    triangles_out = triangles_in;
    model_centroid = mean(node_coords, 1, 'omitnan');

    for triangle_index = 1:size(triangles_out, 1)
        tri = triangles_out(triangle_index, :);
        x = node_coords(tri, :);

        area_vector = cross( ...
            x(2, :) - x(1, :), ...
            x(3, :) - x(1, :));
        triangle_centroid = mean(x, 1);
        outward_hint = triangle_centroid - model_centroid;

        if dot(area_vector, outward_hint) < 0
            triangles_out(triangle_index, [2, 3]) = ...
                triangles_out(triangle_index, [3, 2]);
        end
    end
end

function [node_normals, node_areas] = ...
    compute_pressure_node_geometry( ...
        pressure_node_rows, pressure_triangles, node_coords)

    N_model_nodes = size(node_coords, 1);
    accumulated_area_vector = zeros(N_model_nodes, 3);
    accumulated_area = zeros(N_model_nodes, 1);

    for triangle_index = 1:size(pressure_triangles, 1)
        tri = pressure_triangles(triangle_index, :);
        x = node_coords(tri, :);

        area_vector = 0.5 * cross( ...
            x(2, :) - x(1, :), ...
            x(3, :) - x(1, :));
        area_value = norm(area_vector);

        for local_node = 1:3
            node_row = tri(local_node);
            accumulated_area_vector(node_row, :) = ...
                accumulated_area_vector(node_row, :) + area_vector;
            accumulated_area(node_row) = ...
                accumulated_area(node_row) + area_value / 3;
        end
    end

    node_normals = accumulated_area_vector(pressure_node_rows, :);
    normal_lengths = vecnorm(node_normals, 2, 2);
    valid_normals = normal_lengths > eps;

    node_normals(valid_normals, :) = ...
        node_normals(valid_normals, :) ./ normal_lengths(valid_normals);
    node_normals(~valid_normals, :) = NaN;
    node_areas = accumulated_area(pressure_node_rows);
end

function D = build_edge_difference_matrix( ...
    pressure_triangles, ...
    pressure_node_index_of_model_row, ...
    N_pressure_nodes)

    all_edges_model_rows = [ ...
        pressure_triangles(:, [1, 2]); ...
        pressure_triangles(:, [2, 3]); ...
        pressure_triangles(:, [3, 1])];

    all_edges_model_rows = sort(all_edges_model_rows, 2);
    unique_edges_model_rows = unique(all_edges_model_rows, 'rows');

    edge_node_1 = pressure_node_index_of_model_row( ...
        unique_edges_model_rows(:, 1));
    edge_node_2 = pressure_node_index_of_model_row( ...
        unique_edges_model_rows(:, 2));

    valid_edges = ...
        edge_node_1 > 0 & edge_node_2 > 0 & edge_node_1 ~= edge_node_2;
    edge_node_1 = edge_node_1(valid_edges);
    edge_node_2 = edge_node_2(valid_edges);

    N_edges = numel(edge_node_1);
    D = sparse( ...
        [(1:N_edges)'; (1:N_edges)'], ...
        [edge_node_1; edge_node_2], ...
        [ones(N_edges, 1); -ones(N_edges, 1)], ...
        N_edges, N_pressure_nodes);
end
