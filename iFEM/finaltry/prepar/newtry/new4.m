% =========================================================================
% 完整节点位移 -> 等效节点力 -> 曲面压力
%
% 计算链条：
%
%   1) 从 Abaqus ODB 导出的 CSV 读取每个节点的
%          U1,U2,U3,UR1,UR2,UR3
%
%   2) 使用与 ODB 完全对应的线性有限元刚度矩阵：
%
%          F_equivalent = K_raw * U_full
%
%      对自由节点，F_equivalent 可解释为外部等效节点载荷；
%      对固定节点，F_equivalent 主要是支撑反力。
%
%   3) 在指定的压力加载曲面上建立一致压力载荷矩阵 G：
%
%          F_pressure = G * p
%
%      p 是曲面节点压力。三节点三角形采用线性压力插值：
%
%          f_e = integral(N' * p * n dA)
%
%      对每个三角形，对应的一致载荷系数为：
%
%          A/12 * [2 1 1;
%                  1 2 1;
%                  1 1 2]
%
%   4) 使用自由节点的三个平移力分量反演压力：
%
%          min ||G_fit*p - F_fit||^2
%              + lambda^2 ||D*p||^2
%
%      D 是压力节点相邻边差分矩阵，用于限制压力出现剧烈振荡。
%
% 重要适用条件：
%   A. ODB 与 MTX 必须来自同一模型、同一材料、同一厚度和同一编号；
%   B. 适用于线性静力分析，建议 NLGEOM=OFF；
%   C. ODB 位移必须是全局坐标系；
%   D. CSV 必须包含完整六自由度，不能只提供 U1/U2/U3；
%   E. 反演曲面上最好只存在压力载荷，不应混有集中力、接触力或其他载荷；
%   F. 固定节点的 KU 包含反力，不能直接当作压力节点力参与拟合；
%   G. 当前代码使用三角形曲面。如果原模型是四边形壳单元，
%      read_abaqus_surface_mesh 应将其按一致方式拆成三角形。
%
% 本脚本沿用你现有工程中的：
%
%       read_abaqus_surface_mesh.m
%
% 该函数必须能够返回：
%
%       [node_ids, node_coords, exterior_triangles]
%
% exterior_triangles 中保存 node_coords/node_ids 的行号，而不是 NodeID。
% =========================================================================

clear;
clc;
close all;

%% 1. 用户设置

inp_filepath = 'Zhimian_addnodes_cut.inp';
mtx_filepath = 'Job-5-1_STIF2.mtx';

% 由 export_nodal_u_ur.py 从 ODB 导出的完整节点位移。
% 列顺序必须为：
% NodeID,U1,U2,U3,UR1,UR2,UR3
full_displacement_filepath = 'Abaqus_All_Nodal_U_UR.csv';

% 固定节点文件，一列或多列均可，内容为 INP NodeID。
fixed_filepath = 'fixed_nodes.csv';

% 压力实际施加曲面的节点 NodeID，一列。
% 建议必须提供，以免误把模型全部外表面都当成加载面。
pressure_surface_nodes_filepath = 'pressure_surface_nodes.csv';

% 若压力曲面节点文件不存在：
% false：直接报错，最安全；
% true ：使用模型全部外表面节点，仅适合压力确实作用于整个外表面的情况。
use_all_exterior_if_surface_file_missing = false;

% INP/HyperMesh 坐标 -> Abaqus ODB 全局坐标。
% 只转换节点坐标，不转换从 ODB 直接导出的位移分量。
R_inp_to_odb = diag([1, -1, -1]);

num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];
rotation_dofs = [4, 5, 6];
fixed_dof_components = 1:6;

% 正压力方向：
% -1：沿外法向反方向，即压力压向结构内部；
% +1：沿外法向。
pressure_sign = -1;

% 是否使用模型质心规则，把三角形法向尽量调整为外法向。
% 对封闭或近似凸曲面通常有效；复杂凹曲面需人工检查法向。
orient_surface_normals_outward = true;

% 如果最终压力正负整体颠倒，优先检查三角形法向；
% 也可以把 pressure_sign 从 -1 改成 +1。

% 压力正则化参数选择：
% 'lcurve'：自动使用 L-curve；
% 'fixed' ：使用 fixed_lambda。
lambda_mode = 'lcurve';
fixed_lambda = 1e-3;

% L-curve 搜索范围。这里的 G 已进行尺度归一化，因此 lambda 无量纲。
lambda_candidates = logspace(-8, 2, 90);

% 极小岭稳定项，避免正规方程接近奇异。
pressure_ridge = 1e-12;

% 若完整位移 CSV 某个 MTX 自由度缺失：
% false：报错，推荐；
% true ：缺失自由度按 0 处理，只用于排查文件问题。
allow_missing_displacement_as_zero = false;

displacement_unit = 'mm';
force_unit = 'N';

% 若刚度单位为 N/mm、坐标为 mm，则压力单位为 N/mm^2 = MPa。
pressure_unit = 'MPa';

save_figures = true;

%% 2. 读取节点、坐标和外表面三角形

[node_ids, node_coords, exterior_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath);

node_ids = node_ids(:);
node_coords = node_coords(:,1:3) * R_inp_to_odb;

if isempty(exterior_triangles)
    error('没有读取到外表面三角形。');
end

surface_node_rows = unique(exterior_triangles(:));
surface_node_ids = node_ids(surface_node_rows);
surface_coords = node_coords(surface_node_rows,:);

fprintf('读取模型节点：%d\n', numel(node_ids));
fprintf('模型外表面节点：%d\n', numel(surface_node_rows));
fprintf('模型外表面三角形：%d\n', size(exterior_triangles,1));

%% 3. 读取并组装完整刚度矩阵 K_raw

try
    mtx_data = readmatrix(mtx_filepath, 'FileType', 'text');
catch
    mtx_data = load(mtx_filepath);
end

mtx_data = mtx_data(all(isfinite(mtx_data),2),:);

if size(mtx_data,2) < 5
    error(['MTX 文件至少应包含5列：', ...
           'node_i,dof_i,node_j,dof_j,value。']);
end

node_i = round(mtx_data(:,1));
dof_i  = round(mtx_data(:,2));
node_j = round(mtx_data(:,3));
dof_j  = round(mtx_data(:,4));
values = mtx_data(:,5);

matrix_dof_components = unique([dof_i; dof_j]);

if any(~ismember(matrix_dof_components,1:num_dofs_per_node))
    error('MTX 中出现了超出每节点6自由度范围的自由度编号。');
end

mtx_node_ids = unique([node_i; node_j]);
N_inp_nodes = numel(node_ids);

% node_row_to_mtx_id(r)：
% 第 r 个 INP 节点在 MTX 中使用的节点编号。
node_row_to_mtx_id = nan(N_inp_nodes,1);

if all(ismember(mtx_node_ids,node_ids))
    [~,mapped_rows] = ismember(mtx_node_ids,node_ids);
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = INP NodeID';

elseif all(mtx_node_ids >= 1) && ...
       all(mtx_node_ids <= N_inp_nodes) && ...
       all(mtx_node_ids == round(mtx_node_ids))

    mapped_rows = mtx_node_ids;
    node_row_to_mtx_id(mapped_rows) = mtx_node_ids;
    node_mapping_mode = 'MTX NodeID = INP节点行号';

else
    error(['无法识别 MTX NodeID 与 INP 节点的关系。', newline, ...
           'MTX NodeID 既不是 INP NodeID，也不是有效节点行号。']);
end

row_idx = num_dofs_per_node .* (node_i - 1) + dof_i;
col_idx = num_dofs_per_node .* (node_j - 1) + dof_j;
max_dof = max([row_idx; col_idx]);

K_input = sparse(row_idx,col_idx,values,max_dof,max_dof);

has_upper = nnz(triu(K_input,1)) > 0;
has_lower = nnz(tril(K_input,-1)) > 0;

if xor(has_upper,has_lower)
    % 输入只有上三角或下三角。
    K_raw = K_input + K_input' ...
        - spdiags(diag(K_input),0,max_dof,max_dof);
else
    symmetry_error = norm(K_input-K_input','fro') / ...
        max(norm(K_input,'fro'),eps);

    if symmetry_error > 1e-8
        warning('输入刚度矩阵不是严格对称矩阵，当前进行对称化。');
    end

    K_raw = (K_input + K_input') / 2;
end

[nonzero_row,nonzero_col] = find(K_raw);
existing_dofs = unique([nonzero_row;nonzero_col]);

fprintf('\n========== 刚度矩阵信息 ==========\n');
fprintf('节点映射模式：%s\n',node_mapping_mode);
fprintf('MTX节点数量：%d\n',numel(mtx_node_ids));
fprintf('K_raw尺寸：%d x %d\n',size(K_raw,1),size(K_raw,2));
fprintf('实际存在自由度：%d\n',numel(existing_dofs));
fprintf('非零项数量：%d\n',nnz(K_raw));
fprintf('==================================\n');

%% 4. 读取固定节点并建立固定自由度集合

fixed_data = readmatrix(fixed_filepath);
fixed_nodes_input = unique(round(fixed_data(isfinite(fixed_data))));
fixed_nodes_input = fixed_nodes_input(fixed_nodes_input > 0);
fixed_nodes_input = fixed_nodes_input(:);

[fixed_found,fixed_node_rows] = ...
    ismember(fixed_nodes_input,node_ids);

if any(~fixed_found)
    error('%d 个固定节点不在 INP 节点中。',nnz(~fixed_found));
end

fixed_nodes_mtx = node_row_to_mtx_id(fixed_node_rows);

if any(~isfinite(fixed_nodes_mtx))
    error('部分固定节点没有 MTX 映射。');
end

fixed_global_dofs_matrix = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* (fixed_nodes_mtx(:)-1), ...
    fixed_dof_components);

fixed_global_dofs = fixed_global_dofs_matrix(:);
fixed_existing_dofs = intersect(fixed_global_dofs,existing_dofs);
free_existing_dofs = setdiff(existing_dofs,fixed_existing_dofs);

fprintf('\n========== 固定边界信息 ==========\n');
fprintf('固定节点数量：%d\n',numel(fixed_nodes_input));
fprintf('固定自由度数量：%d\n',numel(fixed_existing_dofs));
fprintf('自由活动自由度：%d\n',numel(free_existing_dofs));
fprintf('==================================\n');

%% 5. 读取完整六自由度位移并组装 U_full

full_displacement_data = readmatrix(full_displacement_filepath);

if size(full_displacement_data,2) < 7
    error(['完整位移文件必须至少有7列：', ...
           'NodeID,U1,U2,U3,UR1,UR2,UR3。']);
end

displacement_ids_raw = full_displacement_data(:,1);
U6_raw = full_displacement_data(:,2:7);

valid_id_rows = isfinite(displacement_ids_raw);

displacement_ids_raw = round( ...
    displacement_ids_raw(valid_id_rows));
U6_raw = U6_raw(valid_id_rows,:);

% 对完全有限的行正常读取。
% 缺少分量的行暂时保留，稍后给出明确诊断。
[displacement_ids_unique,~,displacement_group] = ...
    unique(displacement_ids_raw);

U6_unique = nan(numel(displacement_ids_unique),6);

for component_index = 1:6
    for group_index = 1:numel(displacement_ids_unique)
        current_values = U6_raw( ...
            displacement_group == group_index,component_index);
        current_values = current_values(isfinite(current_values));

        if ~isempty(current_values)
            U6_unique(group_index,component_index) = ...
                mean(current_values);
        end
    end
end

[displacement_found,displacement_node_rows] = ...
    ismember(displacement_ids_unique,node_ids);

if any(~displacement_found)
    warning('%d 个位移 NodeID 不在 INP 中，已忽略。', ...
        nnz(~displacement_found));
end

displacement_ids = ...
    displacement_ids_unique(displacement_found);
displacement_node_rows = ...
    displacement_node_rows(displacement_found);
U6_data = U6_unique(displacement_found,:);

U_full = zeros(max_dof,1);
U_dof_has_value = false(max_dof,1);

for data_index = 1:numel(displacement_ids)

    node_row = displacement_node_rows(data_index);
    mtx_node_id = node_row_to_mtx_id(node_row);

    if ~isfinite(mtx_node_id)
        continue;
    end

    absolute_dofs = ...
        num_dofs_per_node * (mtx_node_id-1) + (1:6);

    for component_index = 1:6
        current_dof = absolute_dofs(component_index);

        if current_dof > max_dof
            continue;
        end

        current_value = U6_data(data_index,component_index);

        if isfinite(current_value)
            U_full(current_dof) = current_value;
            U_dof_has_value(current_dof) = true;
        end
    end
end

missing_existing_dofs = existing_dofs(~U_dof_has_value(existing_dofs));

if ~isempty(missing_existing_dofs)
    fprintf('\n完整位移中缺失的 MTX 自由度数量：%d\n', ...
        numel(missing_existing_dofs));

    print_count = min(20,numel(missing_existing_dofs));
    fprintf('前%d个缺失绝对自由度：',print_count);
    fprintf('%d ',missing_existing_dofs(1:print_count));
    fprintf('\n');

    if allow_missing_displacement_as_zero
        warning('缺失位移自由度将按0处理，KU结果可能不可靠。');
    else
        error(['完整位移文件没有覆盖全部 MTX 自由度。',newline, ...
               '请确认 ODB 已导出 U1,U2,U3,UR1,UR2,UR3。']);
    end
end

fprintf('\n========== 完整位移信息 ==========\n');
fprintf('CSV有效节点数量：%d\n',numel(displacement_ids));
fprintf('覆盖MTX自由度：%d / %d\n', ...
    nnz(U_dof_has_value(existing_dofs)),numel(existing_dofs));
fprintf('完整位移向量范数：%.6e\n',norm(U_full(existing_dofs)));
fprintf('==================================\n');

%% 6. 根据完整位移计算等效节点广义力 F = K*U

F_equivalent_full = K_raw * U_full;

% 提取每个 INP 节点的：
% FX,FY,FZ,MX,MY,MZ。
F_node_6D = nan(N_inp_nodes,6);
U_node_6D = nan(N_inp_nodes,6);

for node_row = 1:N_inp_nodes

    mtx_node_id = node_row_to_mtx_id(node_row);

    if ~isfinite(mtx_node_id)
        continue;
    end

    absolute_dofs = ...
        num_dofs_per_node * (mtx_node_id-1) + (1:6);

    if all(absolute_dofs <= max_dof)
        F_node_6D(node_row,:) = ...
            F_equivalent_full(absolute_dofs)';
        U_node_6D(node_row,:) = ...
            U_full(absolute_dofs)';
    end
end

force_magnitude = vecnorm(F_node_6D(:,1:3),2,2);
moment_magnitude = vecnorm(F_node_6D(:,4:6),2,2);

is_fixed_node = ismember(node_ids,fixed_nodes_input);

% 固定节点上的 KU 主要解释为反力；
% 非固定节点上的 KU 主要解释为外部等效节点力。
free_force_magnitude = force_magnitude;
free_force_magnitude(is_fixed_node) = NaN;

fixed_reaction_sum = sum( ...
    F_node_6D(is_fixed_node,1:3),1,'omitnan');

free_external_force_sum = sum( ...
    F_node_6D(~is_fixed_node,1:3),1,'omitnan');

fprintf('\n========== KU 节点力诊断 ==========\n');
fprintf('固定节点反力合计：[%+.6e, %+.6e, %+.6e] %s\n', ...
    fixed_reaction_sum(1),fixed_reaction_sum(2), ...
    fixed_reaction_sum(3),force_unit);
fprintf('自由节点外力合计：[%+.6e, %+.6e, %+.6e] %s\n', ...
    free_external_force_sum(1),free_external_force_sum(2), ...
    free_external_force_sum(3),force_unit);
fprintf('二者合计：[%+.6e, %+.6e, %+.6e] %s\n', ...
    fixed_reaction_sum(1)+free_external_force_sum(1), ...
    fixed_reaction_sum(2)+free_external_force_sum(2), ...
    fixed_reaction_sum(3)+free_external_force_sum(3), ...
    force_unit);
fprintf('==================================\n');

%% 7. 确定压力加载曲面

if isfile(pressure_surface_nodes_filepath)

    pressure_surface_data = ...
        readmatrix(pressure_surface_nodes_filepath);

    pressure_surface_ids = unique(round( ...
        pressure_surface_data(isfinite(pressure_surface_data))));

    pressure_surface_ids = ...
        pressure_surface_ids(pressure_surface_ids > 0);

    [pressure_id_found,pressure_surface_node_rows_input] = ...
        ismember(pressure_surface_ids,node_ids);

    if any(~pressure_id_found)
        warning('%d 个压力曲面 NodeID 不在 INP 中，已剔除。', ...
            nnz(~pressure_id_found));
    end

    pressure_surface_node_rows_input = ...
        pressure_surface_node_rows_input(pressure_id_found);

else
    if use_all_exterior_if_surface_file_missing
        warning(['没有找到 pressure_surface_nodes.csv，', ...
                 '当前使用全部外表面节点作为压力曲面。']);
        pressure_surface_node_rows_input = surface_node_rows;
    else
        error([ ...
            '没有找到压力曲面节点文件：', ...
            pressure_surface_nodes_filepath,newline, ...
            '请创建一列 NodeID 的 pressure_surface_nodes.csv，', ...
            '或把 use_all_exterior_if_surface_file_missing 设为 true。']);
    end
end

% 只有三个顶点都属于压力曲面节点集的三角形才保留。
is_pressure_triangle = all( ...
    ismember(exterior_triangles, ...
             pressure_surface_node_rows_input),2);

pressure_triangles = exterior_triangles(is_pressure_triangle,:);

if isempty(pressure_triangles)
    error(['压力节点集无法形成三角形曲面。',newline, ...
           '请确认 pressure_surface_nodes.csv 包含完整加载面节点。']);
end

if orient_surface_normals_outward
    pressure_triangles = orient_triangles_by_model_centroid( ...
        pressure_triangles,node_coords);
end

pressure_node_rows = unique(pressure_triangles(:));
pressure_node_ids = node_ids(pressure_node_rows);
pressure_coords = node_coords(pressure_node_rows,:);

N_pressure_nodes = numel(pressure_node_rows);
N_pressure_triangles = size(pressure_triangles,1);

pressure_node_mtx_ids = ...
    node_row_to_mtx_id(pressure_node_rows);

if any(~isfinite(pressure_node_mtx_ids))
    error('压力曲面中存在没有 MTX 映射的节点。');
end

pressure_node_is_fixed = ...
    ismember(pressure_node_rows,fixed_node_rows);

fprintf('\n========== 压力曲面信息 ==========\n');
fprintf('压力曲面节点：%d\n',N_pressure_nodes);
fprintf('压力曲面三角形：%d\n',N_pressure_triangles);
fprintf('压力曲面中的固定节点：%d\n', ...
    nnz(pressure_node_is_fixed));
fprintf('==================================\n');

%% 8. 建立一致压力载荷矩阵 G

% pressure_node_index_of_model_row(r)：
% 模型第 r 个节点在压力未知量中的列号。
pressure_node_index_of_model_row = ...
    zeros(N_inp_nodes,1);

pressure_node_index_of_model_row(pressure_node_rows) = ...
    (1:N_pressure_nodes)';

% 每个三角形：
% 3个受力节点 × 3个压力节点 × 3个力分量 = 27项。
estimated_entries = 27 * N_pressure_triangles;

G_rows = zeros(estimated_entries,1);
G_cols = zeros(estimated_entries,1);
G_vals = zeros(estimated_entries,1);
entry_cursor = 0;

% 三节点线性三角形：
% integral(N_i*N_j dA) / A。
consistent_triangle_weights = ...
    [2,1,1; ...
     1,2,1; ...
     1,1,2] / 12;

triangle_area_vectors = zeros(N_pressure_triangles,3);
triangle_areas = zeros(N_pressure_triangles,1);

for triangle_index = 1:N_pressure_triangles

    tri_rows = pressure_triangles(triangle_index,:);
    tri_coords = node_coords(tri_rows,:);

    area_vector = 0.5 * cross( ...
        tri_coords(2,:)-tri_coords(1,:), ...
        tri_coords(3,:)-tri_coords(1,:));

    triangle_area = norm(area_vector);

    if triangle_area <= eps
        error('压力曲面第%d个三角形面积为0。', ...
            triangle_index);
    end

    triangle_area_vectors(triangle_index,:) = area_vector;
    triangle_areas(triangle_index) = triangle_area;

    local_pressure_indices = ...
        pressure_node_index_of_model_row(tri_rows);

    for local_force_node = 1:3

        force_node_row = tri_rows(local_force_node);
        force_mtx_id = node_row_to_mtx_id(force_node_row);

        force_absolute_dofs = ...
            num_dofs_per_node * (force_mtx_id-1) + ...
            translation_dofs;

        for local_pressure_node = 1:3

            pressure_column = ...
                local_pressure_indices(local_pressure_node);

            coefficient_vector = ...
                pressure_sign * ...
                consistent_triangle_weights( ...
                    local_force_node,local_pressure_node) * ...
                area_vector;

            current_entries = entry_cursor + (1:3);

            G_rows(current_entries) = force_absolute_dofs(:);
            G_cols(current_entries) = pressure_column;
            G_vals(current_entries) = coefficient_vector(:);

            entry_cursor = entry_cursor + 3;
        end
    end
end

G_full = sparse( ...
    G_rows(1:entry_cursor), ...
    G_cols(1:entry_cursor), ...
    G_vals(1:entry_cursor), ...
    max_dof,N_pressure_nodes);

% 节点法向与分摊面积，仅用于输出和显示。
[pressure_node_normals,pressure_node_areas] = ...
    compute_pressure_node_geometry( ...
        pressure_node_rows,pressure_triangles,node_coords);

%% 9. 选择用于压力拟合的节点力分量

pressure_translation_dofs_matrix = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* (pressure_node_mtx_ids(:)-1), ...
    translation_dofs);

% 固定节点 KU 中混入支撑反力，因此不作为压力力拟合行。
fit_pressure_node_mask = ...
    ~pressure_node_is_fixed & ...
    all(ismember(pressure_translation_dofs_matrix, ...
                 existing_dofs),2);

fit_pressure_node_rows = ...
    pressure_node_rows(fit_pressure_node_mask);

fit_pressure_node_ids = ...
    pressure_node_ids(fit_pressure_node_mask);

force_fit_dofs_matrix = ...
    pressure_translation_dofs_matrix(fit_pressure_node_mask,:);

force_fit_dofs = reshape( ...
    force_fit_dofs_matrix',[],1);

G_fit = G_full(force_fit_dofs,:);
F_fit = F_equivalent_full(force_fit_dofs);

% 删除在压力映射中完全为0的方程行。
row_has_pressure_sensitivity = ...
    full(sum(abs(G_fit),2)) > 0;

G_fit = G_fit(row_has_pressure_sensitivity,:);
F_fit = F_fit(row_has_pressure_sensitivity);
force_fit_dofs = force_fit_dofs(row_has_pressure_sensitivity);

if isempty(force_fit_dofs)
    error('没有可用于压力反演的自由节点平移力方程。');
end

% 检查自由节点外力中有多少位于目标压力曲面之外。
all_free_translation_dofs = [];

mapped_nonfixed_rows = find( ...
    isfinite(node_row_to_mtx_id) & ~is_fixed_node);

if ~isempty(mapped_nonfixed_rows)
    mapped_nonfixed_mtx_ids = ...
        node_row_to_mtx_id(mapped_nonfixed_rows);

    all_free_translation_dofs_matrix = bsxfun( ...
        @plus, ...
        num_dofs_per_node .* ...
            (mapped_nonfixed_mtx_ids(:)-1), ...
        translation_dofs);

    all_free_translation_dofs = reshape( ...
        all_free_translation_dofs_matrix',[],1);

    all_free_translation_dofs = intersect( ...
        all_free_translation_dofs,existing_dofs,'stable');
end

non_target_free_translation_dofs = ...
    setdiff(all_free_translation_dofs,force_fit_dofs);

outside_force_norm = norm( ...
    F_equivalent_full(non_target_free_translation_dofs));

all_free_force_norm = norm( ...
    F_equivalent_full(all_free_translation_dofs));

outside_force_ratio = outside_force_norm / ...
    max(all_free_force_norm,eps);

fprintf('\n========== 压力拟合方程 ==========\n');
fprintf('参与拟合的压力曲面自由节点：%d\n', ...
    numel(fit_pressure_node_ids));
fprintf('参与拟合的力分量方程：%d\n',numel(F_fit));
fprintf('压力未知量：%d\n',N_pressure_nodes);
fprintf('目标曲面外自由节点力占比：%.2f %%\n', ...
    100*outside_force_ratio);
fprintf('==================================\n');

if outside_force_ratio > 0.1
    warning([ ...
        '超过10%%的自由节点等效力位于目标压力曲面之外。', ...
        '可能存在其他载荷、压力曲面选择不完整，', ...
        '或 ODB 与 MTX 不完全对应。']);
end

%% 10. 构造压力空间平滑矩阵 D

D_pressure = build_edge_difference_matrix( ...
    pressure_triangles, ...
    pressure_node_index_of_model_row, ...
    N_pressure_nodes);

% 对 G 和 F 同时按 G 的谱范数缩放。
% 这样不改变无正则解，但 lambda 更容易解释。
G_scale = normest(G_fit);

if ~isfinite(G_scale) || G_scale <= eps
    error('压力映射矩阵 G_fit 几乎为零。');
end

G_normalized = G_fit / G_scale;
F_normalized = F_fit / G_scale;

GtG = G_normalized' * G_normalized;
GtF = G_normalized' * F_normalized;
DtD = D_pressure' * D_pressure;

%% 11. 选择正则化参数并求节点压力

switch lower(lambda_mode)

    case 'lcurve'

        N_lambda = numel(lambda_candidates);
        residual_norms = zeros(N_lambda,1);
        regularization_norms = zeros(N_lambda,1);
        pressure_solutions = cell(N_lambda,1);

        for lambda_index = 1:N_lambda

            current_lambda = ...
                lambda_candidates(lambda_index);

            pressure_system = ...
                GtG + ...
                current_lambda^2 * DtD + ...
                pressure_ridge * speye(N_pressure_nodes);

            current_pressure = pressure_system \ GtF;

            pressure_solutions{lambda_index} = ...
                current_pressure;

            residual_norms(lambda_index) = norm( ...
                G_normalized*current_pressure - ...
                F_normalized);

            regularization_norms(lambda_index) = sqrt( ...
                norm(D_pressure*current_pressure)^2 + ...
                pressure_ridge*norm(current_pressure)^2);
        end

        log_lambda = log10(lambda_candidates(:));
        log_residual = log10(max(residual_norms,eps));
        log_regularization = ...
            log10(max(regularization_norms,eps));

        dx = gradient(log_residual,log_lambda);
        dy = gradient(log_regularization,log_lambda);
        d2x = gradient(dx,log_lambda);
        d2y = gradient(dy,log_lambda);

        curvature = abs( ...
            (dx.*d2y - dy.*d2x) ./ ...
            max((dx.^2 + dy.^2).^(3/2),eps));

        curvature(~isfinite(curvature)) = -inf;

        valid_lambda_indices = 5:(N_lambda-4);

        if isempty(valid_lambda_indices)
            valid_lambda_indices = 1:N_lambda;
        end

        [~,best_offset] = ...
            max(curvature(valid_lambda_indices));

        best_lambda_index = ...
            valid_lambda_indices(best_offset);

        lambda_selected = ...
            lambda_candidates(best_lambda_index);

        pressure_nodal = ...
            pressure_solutions{best_lambda_index};

        figure( ...
            'Name','Pressure L-curve', ...
            'Color','w', ...
            'Position',[100,100,950,420]);

        subplot(1,2,1);
        loglog(residual_norms,regularization_norms, ...
            'LineWidth',1.5);
        hold on;
        loglog( ...
            residual_norms(best_lambda_index), ...
            regularization_norms(best_lambda_index), ...
            'ro','MarkerFaceColor','r');
        grid on;
        xlabel('||G p - F||');
        ylabel('||D p||');
        title('压力反演 L-curve');

        subplot(1,2,2);
        semilogx(lambda_candidates,curvature, ...
            'LineWidth',1.5);
        hold on;
        semilogx(lambda_selected, ...
            curvature(best_lambda_index), ...
            'ro','MarkerFaceColor','r');
        grid on;
        xlabel('\lambda');
        ylabel('L-curve 曲率');
        title(sprintf('选取 \\lambda = %.3e', ...
            lambda_selected));

        if save_figures
            exportgraphics(gcf, ...
                'Pressure_Inversion_Lcurve.png', ...
                'Resolution',300);
        end

    case 'fixed'

        lambda_selected = fixed_lambda;

        pressure_system = ...
            GtG + ...
            lambda_selected^2 * DtD + ...
            pressure_ridge * speye(N_pressure_nodes);

        pressure_nodal = pressure_system \ GtF;

    otherwise
        error('lambda_mode 必须是 lcurve 或 fixed。');
end

%% 12. 根据压力重构节点力并计算残差

F_pressure_full = G_full * pressure_nodal;

F_pressure_fit = ...
    F_pressure_full(force_fit_dofs);

F_equivalent_fit = ...
    F_equivalent_full(force_fit_dofs);

force_fit_residual = ...
    F_pressure_fit - F_equivalent_fit;

relative_force_fit_error = ...
    norm(force_fit_residual) / ...
    max(norm(F_equivalent_fit),eps);

force_fit_cosine = dot( ...
    F_pressure_fit,F_equivalent_fit) / ...
    max(norm(F_pressure_fit)* ...
        norm(F_equivalent_fit),eps);

pressure_force_node_3D = ...
    nan(N_pressure_nodes,3);

KU_force_pressure_node_3D = ...
    nan(N_pressure_nodes,3);

for pressure_index = 1:N_pressure_nodes

    current_mtx_id = ...
        pressure_node_mtx_ids(pressure_index);

    current_translation_dofs = ...
        num_dofs_per_node * (current_mtx_id-1) + ...
        translation_dofs;

    pressure_force_node_3D(pressure_index,:) = ...
        F_pressure_full(current_translation_dofs)';

    KU_force_pressure_node_3D(pressure_index,:) = ...
        F_equivalent_full(current_translation_dofs)';
end

pressure_force_magnitude = ...
    vecnorm(pressure_force_node_3D,2,2);

KU_force_pressure_magnitude = ...
    vecnorm(KU_force_pressure_node_3D,2,2);

pressure_total_force = ...
    sum(pressure_force_node_3D,1,'omitnan');

fprintf('\n========== 压力反演结果 ==========\n');
fprintf('选取 lambda：%.6e\n',lambda_selected);
fprintf('节点压力最小值：%.6e %s\n', ...
    min(pressure_nodal),pressure_unit);
fprintf('节点压力最大值：%.6e %s\n', ...
    max(pressure_nodal),pressure_unit);
fprintf('节点压力平均值：%.6e %s\n', ...
    mean(pressure_nodal),pressure_unit);
fprintf('压力力拟合相对误差：%.2f %%\n', ...
    100*relative_force_fit_error);
fprintf('压力力与KU力余弦相似度：%.6f\n', ...
    force_fit_cosine);
fprintf('重构压力总力：[%+.6e, %+.6e, %+.6e] %s\n', ...
    pressure_total_force(1), ...
    pressure_total_force(2), ...
    pressure_total_force(3),force_unit);
fprintf('==================================\n');

if mean(pressure_nodal) < 0
    warning([ ...
        '平均压力为负。若物理上压力应为正，', ...
        '请检查曲面三角形法向和 pressure_sign。']);
end

%% 13. 输出完整节点 KU 力表

is_pressure_node = ...
    ismember(node_ids,pressure_node_ids);

equivalent_force_table = table( ...
    node_ids(:), ...
    node_coords(:,1), ...
    node_coords(:,2), ...
    node_coords(:,3), ...
    U_node_6D(:,1), ...
    U_node_6D(:,2), ...
    U_node_6D(:,3), ...
    U_node_6D(:,4), ...
    U_node_6D(:,5), ...
    U_node_6D(:,6), ...
    F_node_6D(:,1), ...
    F_node_6D(:,2), ...
    F_node_6D(:,3), ...
    F_node_6D(:,4), ...
    F_node_6D(:,5), ...
    F_node_6D(:,6), ...
    force_magnitude, ...
    moment_magnitude, ...
    is_fixed_node(:), ...
    is_pressure_node(:), ...
    'VariableNames',{ ...
        'NodeID','X','Y','Z', ...
        'U1','U2','U3','UR1','UR2','UR3', ...
        'FX_from_KU','FY_from_KU','FZ_from_KU', ...
        'MX_from_KU','MY_from_KU','MZ_from_KU', ...
        'ForceMagnitude','MomentMagnitude', ...
        'IsFixedNode','IsPressureSurfaceNode'});

writetable( ...
    equivalent_force_table, ...
    'Equivalent_Nodal_Force_From_KU.csv');

fprintf('已保存：Equivalent_Nodal_Force_From_KU.csv\n');

%% 14. 输出压力节点结果表

pressure_result_table = table( ...
    pressure_node_ids(:), ...
    pressure_coords(:,1), ...
    pressure_coords(:,2), ...
    pressure_coords(:,3), ...
    pressure_node_normals(:,1), ...
    pressure_node_normals(:,2), ...
    pressure_node_normals(:,3), ...
    pressure_node_areas(:), ...
    pressure_nodal(:), ...
    KU_force_pressure_node_3D(:,1), ...
    KU_force_pressure_node_3D(:,2), ...
    KU_force_pressure_node_3D(:,3), ...
    KU_force_pressure_magnitude(:), ...
    pressure_force_node_3D(:,1), ...
    pressure_force_node_3D(:,2), ...
    pressure_force_node_3D(:,3), ...
    pressure_force_magnitude(:), ...
    pressure_node_is_fixed(:), ...
    'VariableNames',{ ...
        'NodeID','X','Y','Z', ...
        'NormalX','NormalY','NormalZ', ...
        'NodalArea','Pressure_Reconstructed', ...
        'FX_From_KU','FY_From_KU','FZ_From_KU', ...
        'KU_ForceMagnitude', ...
        'FX_From_Pressure','FY_From_Pressure', ...
        'FZ_From_Pressure','PressureForceMagnitude', ...
        'IsFixedNode'});

writetable( ...
    pressure_result_table, ...
    'Pressure_From_Complete_Displacement.csv');

fprintf('已保存：Pressure_From_Complete_Displacement.csv\n');

%% 15. 输出拟合自由度逐分量比较表

fit_dof_node_mtx_id = ...
    floor((force_fit_dofs-1)/num_dofs_per_node) + 1;

fit_dof_component = ...
    mod(force_fit_dofs-1,num_dofs_per_node) + 1;

% MTX NodeID -> INP NodeID。
fit_dof_inp_node_id = nan(numel(force_fit_dofs),1);

for row_index = 1:numel(force_fit_dofs)
    current_mtx_id = fit_dof_node_mtx_id(row_index);

    matching_model_row = find( ...
        node_row_to_mtx_id == current_mtx_id,1);

    if ~isempty(matching_model_row)
        fit_dof_inp_node_id(row_index) = ...
            node_ids(matching_model_row);
    end
end

force_fit_table = table( ...
    fit_dof_inp_node_id, ...
    fit_dof_component, ...
    F_equivalent_fit, ...
    F_pressure_fit, ...
    force_fit_residual, ...
    'VariableNames',{ ...
        'NodeID','Component', ...
        'Force_From_KU', ...
        'Force_From_ReconstructedPressure', ...
        'Force_Error'});

writetable( ...
    force_fit_table, ...
    'Pressure_Force_Fit_Details.csv');

fprintf('已保存：Pressure_Force_Fit_Details.csv\n');

%% 16. 绘制自由节点 KU 等效力模长云图

node_row_to_surface_row = ...
    zeros(N_inp_nodes,1);

node_row_to_surface_row(surface_node_rows) = ...
    (1:numel(surface_node_rows))';

surface_triangles_local = ...
    node_row_to_surface_row(exterior_triangles);

surface_triangles_local = ...
    surface_triangles_local( ...
        all(surface_triangles_local > 0,2),:);

surface_free_force_magnitude = ...
    free_force_magnitude(surface_node_rows);

figure( ...
    'Name','自由节点 KU 等效力模长', ...
    'Color','w', ...
    'Position',[100,100,900,700]);

trisurf( ...
    surface_triangles_local, ...
    surface_coords(:,1), ...
    surface_coords(:,2), ...
    surface_coords(:,3), ...
    surface_free_force_magnitude, ...
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
title('自由节点等效外力模长 ||F|| = ||KU||');
cb = colorbar;
ylabel(cb,['ForceMagnitude (',force_unit,')']);
colormap(jet(256));

if save_figures
    exportgraphics(gcf, ...
        'Equivalent_Free_Nodal_Force_Cloud.png', ...
        'Resolution',300);
end

%% 17. 绘制重构压力云图

pressure_row_to_local = zeros(N_inp_nodes,1);
pressure_row_to_local(pressure_node_rows) = ...
    (1:N_pressure_nodes)';

pressure_triangles_local = ...
    pressure_row_to_local(pressure_triangles);

figure( ...
    'Name','完整位移反演曲面压力', ...
    'Color','w', ...
    'Position',[100,100,900,700]);

trisurf( ...
    pressure_triangles_local, ...
    pressure_coords(:,1), ...
    pressure_coords(:,2), ...
    pressure_coords(:,3), ...
    pressure_nodal, ...
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
title('由完整六自由度位移重构的曲面压力');
cb = colorbar;
ylabel(cb,['Pressure (',pressure_unit,')']);
colormap(jet(256));

hold on;
scatter3( ...
    pressure_coords(pressure_node_is_fixed,1), ...
    pressure_coords(pressure_node_is_fixed,2), ...
    pressure_coords(pressure_node_is_fixed,3), ...
    25,'k','filled');
hold off;

if save_figures
    exportgraphics(gcf, ...
        'Pressure_From_Complete_Displacement_Cloud.png', ...
        'Resolution',300);
end

%% 18. 显示 KU 等效外力最大的自由节点

free_node_rows = find( ...
    ~is_fixed_node & isfinite(force_magnitude));

[~,force_order] = sort( ...
    force_magnitude(free_node_rows),'descend');

top_count = min(15,numel(force_order));
top_force_rows = ...
    free_node_rows(force_order(1:top_count));

top_force_table = equivalent_force_table( ...
    top_force_rows, ...
    {'NodeID','FX_from_KU','FY_from_KU','FZ_from_KU', ...
     'ForceMagnitude','IsPressureSurfaceNode'});

fprintf('\nKU 等效外力最大的自由节点：\n');
disp(top_force_table);

fprintf('\n程序运行完成。\n');

%% =========================================================================
% 局部函数
% =========================================================================

function triangles_out = orient_triangles_by_model_centroid( ...
    triangles_in,node_coords)

    triangles_out = triangles_in;
    model_centroid = mean(node_coords,1,'omitnan');

    for triangle_index = 1:size(triangles_out,1)

        tri = triangles_out(triangle_index,:);
        x = node_coords(tri,:);

        area_vector = cross( ...
            x(2,:)-x(1,:), ...
            x(3,:)-x(1,:));

        triangle_centroid = mean(x,1);
        outward_hint = triangle_centroid-model_centroid;

        if dot(area_vector,outward_hint) < 0
            triangles_out(triangle_index,[2,3]) = ...
                triangles_out(triangle_index,[3,2]);
        end
    end
end

function [node_normals,node_areas] = ...
    compute_pressure_node_geometry( ...
        pressure_node_rows,pressure_triangles,node_coords)

    N_model_nodes = size(node_coords,1);
    accumulated_area_vector = zeros(N_model_nodes,3);
    accumulated_area = zeros(N_model_nodes,1);

    for triangle_index = 1:size(pressure_triangles,1)

        tri = pressure_triangles(triangle_index,:);
        x = node_coords(tri,:);

        area_vector = 0.5 * cross( ...
            x(2,:)-x(1,:), ...
            x(3,:)-x(1,:));

        area_value = norm(area_vector);

        for local_node = 1:3
            node_row = tri(local_node);

            accumulated_area_vector(node_row,:) = ...
                accumulated_area_vector(node_row,:) + ...
                area_vector;

            accumulated_area(node_row) = ...
                accumulated_area(node_row) + ...
                area_value/3;
        end
    end

    node_normals = ...
        accumulated_area_vector(pressure_node_rows,:);

    normal_lengths = vecnorm(node_normals,2,2);

    valid_normals = normal_lengths > eps;

    node_normals(valid_normals,:) = ...
        node_normals(valid_normals,:) ./ ...
        normal_lengths(valid_normals);

    node_normals(~valid_normals,:) = NaN;

    node_areas = ...
        accumulated_area(pressure_node_rows);
end

function D = build_edge_difference_matrix( ...
    pressure_triangles, ...
    pressure_node_index_of_model_row, ...
    N_pressure_nodes)

    all_edges_model_rows = [ ...
        pressure_triangles(:,[1,2]); ...
        pressure_triangles(:,[2,3]); ...
        pressure_triangles(:,[3,1])];

    all_edges_model_rows = sort( ...
        all_edges_model_rows,2);

    unique_edges_model_rows = unique( ...
        all_edges_model_rows,'rows');

    edge_node_1 = pressure_node_index_of_model_row( ...
        unique_edges_model_rows(:,1));

    edge_node_2 = pressure_node_index_of_model_row( ...
        unique_edges_model_rows(:,2));

    valid_edges = ...
        edge_node_1 > 0 & edge_node_2 > 0 & ...
        edge_node_1 ~= edge_node_2;

    edge_node_1 = edge_node_1(valid_edges);
    edge_node_2 = edge_node_2(valid_edges);

    N_edges = numel(edge_node_1);

    D = sparse( ...
        [(1:N_edges)';(1:N_edges)'], ...
        [edge_node_1;edge_node_2], ...
        [ones(N_edges,1);-ones(N_edges,1)], ...
        N_edges,N_pressure_nodes);
end
