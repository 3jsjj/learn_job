% =========================================================================
% 完整六自由度位移 -> 等效节点力 -> 自动识别曲面压力（无 pressure_nodes）
%
% 计算链条：
%
%   1) 从 Abaqus ODB 导出的 CSV 读取全部节点：
%
%          NodeID,U1,U2,U3,UR1,UR2,UR3
%
%   2) 使用与 ODB 完全对应的线性有限元刚度矩阵：
%
%          F_equivalent = K_raw * U_full
%
%   3) 把“全部外表面节点”都作为候选压力节点，并建立一致压力载荷矩阵：
%
%          F_pressure = G * p
%
%   4) 不提供 pressure_nodes，不预先指定受压区域，而是求解：
%
%          min_{p >= 0}
%              1/2 ||G p - F||^2
%            + 1/2 alpha^2 ||D p||^2
%            + beta ||p||_1
%            + 1/2 gamma ||p||^2
%
%      其中：
%        p >= 0       保证压力不反向；
%        ||D p||^2    保证相邻节点压力不过度振荡；
%        ||p||_1      促使未受压区域自动变为零，自动识别加载区域；
%        gamma        是极小数值稳定项。
%
%   5) 程序同时尝试压力沿外法向和沿外法向反方向，
%      自动选取节点力拟合残差较小的方向。
%
% 输出：
%   Equivalent_Nodal_Force_From_KU.csv
%   Pressure_From_Complete_Displacement_AutoSurface.csv
%   Pressure_Force_Fit_Details_AutoSurface.csv
%   Equivalent_Free_Nodal_Force_Cloud.png
%   Pressure_From_Complete_Displacement_AutoSurface.png
%
% 重要适用条件：
%   A. ODB 与 MTX 必须来自同一个线性模型；
%   B. 建议 NLGEOM=OFF；
%   C. ODB 位移和转角必须为全局坐标系；
%   D. CSV 必须包含完整六自由度；
%   E. 目标外载荷主要应为曲面法向压力；
%   F. 固定节点的 KU 主要是支撑反力，因此不参与压力拟合；
%   G. 若模型同时存在集中力、接触力、重力等，压力拟合残差会升高；
%   H. 本方法可以自动识别受压区域，但这仍是正则化逆问题，
%      稀疏权重和光滑权重会影响加载区域大小与压力平滑程度。
%
% 本脚本继续调用：
%
%       read_abaqus_surface_mesh.m
%
% 函数应返回：
%
%       [node_ids, node_coords, exterior_triangles]
%
% exterior_triangles 保存节点行号，而不是 NodeID。
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

% INP/HyperMesh 坐标 -> Abaqus ODB 全局坐标。
% 只转换节点坐标，不转换从 ODB 直接导出的位移和转角。
R_inp_to_odb = diag([1, -1, -1]);

num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];
rotation_dofs = [4, 5, 6];
fixed_dof_components = 1:6;

% 自动比较两个全局压力方向：
%   -1：沿外法向反方向；
%   +1：沿外法向。
auto_choose_pressure_sign = true;
manual_pressure_sign = -1;
pressure_sign_candidates = [-1, +1];

% 使用模型质心规则，把三角形法向尽量调整为外法向。
% 对复杂凹曲面或非封闭曲面，应结合输出法向进行检查。
orient_surface_normals_outward = true;

% -------------------------------------------------------------------------
% 自动识别受压区域的正则化参数
% -------------------------------------------------------------------------
%
% smoothness_weight：
%   越大，压力越平滑，但峰值会被压低、加载区域会变宽。
%
% sparsity_fraction：
%   beta = sparsity_fraction * beta_max。
%   越大，非受压节点越容易被压成0，加载区域越集中；
%   太大可能漏掉真实低压力区域。
%
% 由于 G 和 D 都会自动归一化，这两个参数是无量纲的。
smoothness_weight = 2e-2;
sparsity_fraction = 5e-3;
ridge_weight = 1e-10;

% FISTA 迭代设置，不依赖 Optimization Toolbox。
pressure_solver_max_iterations = 6000;
pressure_solver_tolerance = 1e-9;

% 后处理：压力达到最大值多少比例时，标记为“活动受压节点”。
% 这只影响活动区域标记，不影响反演求解。
active_pressure_relative_threshold = 0.02;

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

%% 7. 全部外表面自动作为候选压力曲面

% 不读取 pressure_nodes。
% 所有外表面节点和三角形都进入压力反演候选集合。
pressure_triangles = exterior_triangles;

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
    error('部分外表面节点没有 MTX 映射。');
end

pressure_node_is_fixed = ...
    ismember(pressure_node_rows,fixed_node_rows);

fprintf('\n========== 自动候选压力曲面 ==========\n');
fprintf('候选节点：全部外表面节点，共 %d 个\n', ...
    N_pressure_nodes);
fprintf('候选三角形：全部外表面三角形，共 %d 个\n', ...
    N_pressure_triangles);
fprintf('候选曲面中的固定节点：%d\n', ...
    nnz(pressure_node_is_fixed));
fprintf('不使用 pressure_nodes 文件。\n');
fprintf('======================================\n');

%% 8. 在全部外表面建立一致压力载荷矩阵 G

% pressure_node_index_of_model_row(r)：
% 模型第 r 个节点在压力未知量 p 中的列号。
pressure_node_index_of_model_row = zeros(N_inp_nodes,1);
pressure_node_index_of_model_row(pressure_node_rows) = ...
    (1:N_pressure_nodes)';

% 每个三角形：
% 3个受力节点 × 3个压力节点 × 3个力分量 = 27项。
estimated_entries = 27 * N_pressure_triangles;

G_rows = zeros(estimated_entries,1);
G_cols = zeros(estimated_entries,1);
G_vals = zeros(estimated_entries,1);
entry_cursor = 0;

% 三节点线性三角形的一致压力载荷权重：
%
% integral_A N_i N_j dA
%       = A/12 * [2 1 1; 1 2 1; 1 1 2]
consistent_triangle_weights = ...
    [2,1,1; ...
     1,2,1; ...
     1,1,2] / 12;

triangle_area_vectors = zeros(N_pressure_triangles,3);
triangle_areas = zeros(N_pressure_triangles,1);

for triangle_index = 1:N_pressure_triangles

    tri_rows = pressure_triangles(triangle_index,:);
    tri_coords = node_coords(tri_rows,:);

    % area_vector 的模长等于三角形面积，方向由节点顺序确定。
    area_vector = 0.5 * cross( ...
        tri_coords(2,:)-tri_coords(1,:), ...
        tri_coords(3,:)-tri_coords(1,:));

    triangle_area = norm(area_vector);

    if triangle_area <= eps
        error('外表面第%d个三角形面积为0。', ...
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

            % 这里暂不乘正负号，先构造沿外法向的 G。
            coefficient_vector = ...
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

G_outward_full = sparse( ...
    G_rows(1:entry_cursor), ...
    G_cols(1:entry_cursor), ...
    G_vals(1:entry_cursor), ...
    max_dof,N_pressure_nodes);

% 节点法向与分摊面积，用于输出和显示。
[pressure_node_normals,pressure_node_areas] = ...
    compute_pressure_node_geometry( ...
        pressure_node_rows,pressure_triangles,node_coords);

%% 9. 选择全部外表面自由节点的平移力作为拟合方程

pressure_translation_dofs_matrix = bsxfun( ...
    @plus, ...
    num_dofs_per_node .* (pressure_node_mtx_ids(:)-1), ...
    translation_dofs);

% 固定节点的 KU 中包含支撑反力，不作为外部压力力拟合方程。
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

G_outward_fit = G_outward_full(force_fit_dofs,:);
F_fit = F_equivalent_full(force_fit_dofs);

% 删除压力映射完全为0的力方程。
row_has_pressure_sensitivity = ...
    full(sum(abs(G_outward_fit),2)) > 0;

G_outward_fit = ...
    G_outward_fit(row_has_pressure_sensitivity,:);

F_fit = F_fit(row_has_pressure_sensitivity);
force_fit_dofs = force_fit_dofs(row_has_pressure_sensitivity);

if isempty(force_fit_dofs)
    error('没有可用于压力反演的外表面自由节点平移力方程。');
end

% 检查自由节点 KU 力中有多少位于外表面以外。
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

non_surface_free_translation_dofs = ...
    setdiff(all_free_translation_dofs,force_fit_dofs);

non_surface_force_norm = norm( ...
    F_equivalent_full(non_surface_free_translation_dofs));

all_free_force_norm = norm( ...
    F_equivalent_full(all_free_translation_dofs));

non_surface_force_ratio = non_surface_force_norm / ...
    max(all_free_force_norm,eps);

fprintf('\n========== 自动压力拟合方程 ==========\n');
fprintf('参与拟合的外表面自由节点：%d\n', ...
    numel(fit_pressure_node_ids));
fprintf('参与拟合的力分量方程：%d\n',numel(F_fit));
fprintf('候选压力未知量：%d\n',N_pressure_nodes);
fprintf('外表面以外自由节点力占比：%.2f %%\n', ...
    100*non_surface_force_ratio);
fprintf('======================================\n');

if non_surface_force_ratio > 0.1
    warning([ ...
        '超过10%%的自由节点等效力不在外表面拟合方程中。', ...
        '可能存在内部集中力、接触、惯性载荷，', ...
        '或 ODB 与 MTX 不完全对应。']);
end

%% 10. 构造压力光滑矩阵并进行尺度归一化

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

%% 11. 非负稀疏压力反演，并自动选择压力方向

if auto_choose_pressure_sign
    tested_pressure_signs = pressure_sign_candidates(:);
else
    tested_pressure_signs = manual_pressure_sign;
end

N_signs = numel(tested_pressure_signs);

pressure_candidates = cell(N_signs,1);
solver_info_candidates = cell(N_signs,1);
sign_relative_residuals = inf(N_signs,1);
sign_active_counts = zeros(N_signs,1);
sign_l1_weights = zeros(N_signs,1);

fprintf('\n========== 压力方向自动比较 ==========\n');

for sign_index = 1:N_signs

    current_sign = tested_pressure_signs(sign_index);

    G_candidate = ...
        current_sign * G_outward_normalized;

    % 对非负 L1 问题，全部压力变为0的临界 beta 为：
    %
    % beta_max = max(G'F)
    %
    % 取其一定比例，使求解既允许非零压力，又倾向局部加载区域。
    beta_max = max(full(G_candidate' * F_normalized));
    beta_max = max(beta_max,0);

    current_l1_weight = ...
        sparsity_fraction * beta_max;

    [current_pressure,current_info] = ...
        solve_sparse_nonnegative_pressure_fista( ...
            G_candidate, ...
            F_normalized, ...
            D_normalized, ...
            smoothness_weight, ...
            current_l1_weight, ...
            ridge_weight, ...
            pressure_solver_max_iterations, ...
            pressure_solver_tolerance);

    pressure_candidates{sign_index} = current_pressure;
    solver_info_candidates{sign_index} = current_info;
    sign_relative_residuals(sign_index) = ...
        current_info.relative_residual;
    sign_l1_weights(sign_index) = current_l1_weight;

    if max(current_pressure) > 0
        sign_active_counts(sign_index) = nnz( ...
            current_pressure >= ...
            active_pressure_relative_threshold * ...
            max(current_pressure));
    end

    fprintf(['方向 %+d：残差 %.4f，活动节点 %d，', ...
             'beta %.3e，迭代 %d\n'], ...
        current_sign, ...
        sign_relative_residuals(sign_index), ...
        sign_active_counts(sign_index), ...
        current_l1_weight, ...
        current_info.iterations);
end

fprintf('======================================\n');

% 方向选择以力拟合相对残差为主。
[~,best_sign_index] = min(sign_relative_residuals);

pressure_sign_selected = ...
    tested_pressure_signs(best_sign_index);

pressure_nodal = ...
    pressure_candidates{best_sign_index};

pressure_solver_info = ...
    solver_info_candidates{best_sign_index};

sparsity_weight_selected = ...
    sign_l1_weights(best_sign_index);

G_full = pressure_sign_selected * G_outward_full;

if ~any(pressure_nodal > 0)
    warning([ ...
        '自动反演得到的压力全部为0。', ...
        '请检查完整位移、刚度矩阵、坐标变换和法向方向；', ...
        '也可适当减小 sparsity_fraction。']);
end

%% 12. 根据自动压力重构节点力并计算残差

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

if max(pressure_nodal) > 0
    pressure_relative_to_max = ...
        pressure_nodal / max(pressure_nodal);

    is_active_pressure_node = ...
        pressure_relative_to_max >= ...
        active_pressure_relative_threshold;
else
    pressure_relative_to_max = ...
        zeros(size(pressure_nodal));

    is_active_pressure_node = ...
        false(size(pressure_nodal));
end

fprintf('\n========== 自动压力反演结果 ==========\n');
fprintf('自动选取压力方向：%+d\n', ...
    pressure_sign_selected);
fprintf('光滑权重 alpha：%.6e\n', ...
    smoothness_weight);
fprintf('稀疏权重 beta：%.6e\n', ...
    sparsity_weight_selected);
fprintf('FISTA迭代次数：%d\n', ...
    pressure_solver_info.iterations);
fprintf('候选外表面节点：%d\n', ...
    N_pressure_nodes);
fprintf('自动识别活动受压节点：%d\n', ...
    nnz(is_active_pressure_node));
fprintf('活动阈值：最大压力的 %.2f %%\n', ...
    100*active_pressure_relative_threshold);
fprintf('节点压力最小值：%.6e %s\n', ...
    min(pressure_nodal),pressure_unit);
fprintf('节点压力最大值：%.6e %s\n', ...
    max(pressure_nodal),pressure_unit);

if any(is_active_pressure_node)
    active_pressure_mean = ...
        mean(pressure_nodal(is_active_pressure_node),'omitnan');
else
    active_pressure_mean = NaN;
end

fprintf('活动节点平均压力：%.6e %s\n', ...
    active_pressure_mean,pressure_unit);
fprintf('压力力拟合相对误差：%.2f %%\n', ...
    100*relative_force_fit_error);
fprintf('压力力与KU力余弦相似度：%.6f\n', ...
    force_fit_cosine);
fprintf('重构压力总力：[%+.6e, %+.6e, %+.6e] %s\n', ...
    pressure_total_force(1), ...
    pressure_total_force(2), ...
    pressure_total_force(3),force_unit);
fprintf('======================================\n');

%% 13. 输出完整节点 KU 力表

is_exterior_surface_node = ...
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
    is_exterior_surface_node(:), ...
    'VariableNames',{ ...
        'NodeID','X','Y','Z', ...
        'U1','U2','U3','UR1','UR2','UR3', ...
        'FX_from_KU','FY_from_KU','FZ_from_KU', ...
        'MX_from_KU','MY_from_KU','MZ_from_KU', ...
        'ForceMagnitude','MomentMagnitude', ...
        'IsFixedNode','IsExteriorSurfaceNode'});

writetable( ...
    equivalent_force_table, ...
    'Equivalent_Nodal_Force_From_KU.csv');

fprintf('已保存：Equivalent_Nodal_Force_From_KU.csv\n');

%% 14. 输出全部外表面的自动压力结果

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
    pressure_relative_to_max(:), ...
    is_active_pressure_node(:), ...
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
        'Pressure_RelativeToMax', ...
        'IsActivePressureNode', ...
        'FX_From_KU','FY_From_KU','FZ_From_KU', ...
        'KU_ForceMagnitude', ...
        'FX_From_Pressure','FY_From_Pressure', ...
        'FZ_From_Pressure','PressureForceMagnitude', ...
        'IsFixedNode'});

writetable( ...
    pressure_result_table, ...
    'Pressure_From_Complete_Displacement_AutoSurface.csv');

fprintf(['已保存：', ...
    'Pressure_From_Complete_Displacement_AutoSurface.csv\n']);

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
    'Pressure_Force_Fit_Details_AutoSurface.csv');

fprintf('已保存：Pressure_Force_Fit_Details_AutoSurface.csv\n');

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
    'Name','自动识别外表面压力', ...
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
title('由完整六自由度位移自动识别的外表面压力');
cb = colorbar;
ylabel(cb,['Pressure (',pressure_unit,')']);
colormap(jet(256));

hold on;
scatter3( ...
    pressure_coords(pressure_node_is_fixed,1), ...
    pressure_coords(pressure_node_is_fixed,2), ...
    pressure_coords(pressure_node_is_fixed,3), ...
    25,'k','filled');

scatter3( ...
    pressure_coords(is_active_pressure_node,1), ...
    pressure_coords(is_active_pressure_node,2), ...
    pressure_coords(is_active_pressure_node,3), ...
    18,'ko');
hold off;

if save_figures
    exportgraphics(gcf, ...
        'Pressure_From_Complete_Displacement_AutoSurface.png', ...
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
     'ForceMagnitude','IsExteriorSurfaceNode'});

fprintf('\nKU 等效外力最大的自由节点：\n');
disp(top_force_table);

fprintf('\n程序运行完成。\n');

%% =========================================================================
% 局部函数
% =========================================================================


function [pressure,info] = ...
    solve_sparse_nonnegative_pressure_fista( ...
        G,F,D,lambda_smooth,lambda_l1,lambda_ridge, ...
        max_iterations,tolerance)
% -------------------------------------------------------------------------
% 求解：
%
%   min_{p>=0}
%       1/2||Gp-F||^2
%     + 1/2 lambda_smooth^2 ||Dp||^2
%     + lambda_l1 ||p||_1
%     + 1/2 lambda_ridge ||p||^2
%
% 使用 FISTA：
%   光滑部分做梯度下降；
%   非负 L1 项的近端算子为 max(z-lambda_l1/L,0)。
% -------------------------------------------------------------------------

    N_pressure = size(G,2);

    pressure = zeros(N_pressure,1);
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

    lipschitz_constant = ...
        max(lipschitz_constant,eps);

    objective_previous = inf;
    converged = false;

    for iteration = 1:max_iterations

        gradient_value = ...
            G' * (G*extrapolated_pressure-F) + ...
            lambda_smooth^2 * ...
                (D'*(D*extrapolated_pressure)) + ...
            lambda_ridge * extrapolated_pressure;

        proximal_argument = ...
            extrapolated_pressure - ...
            gradient_value/lipschitz_constant;

        pressure_new = max( ...
            proximal_argument - ...
            lambda_l1/lipschitz_constant, ...
            0);

        momentum_new = ...
            (1+sqrt(1+4*momentum^2))/2;

        extrapolated_new = ...
            pressure_new + ...
            ((momentum-1)/momentum_new) * ...
            (pressure_new-pressure);

        current_residual = G*pressure_new-F;

        current_objective = ...
            0.5*norm(current_residual)^2 + ...
            0.5*lambda_smooth^2 * ...
                norm(D*pressure_new)^2 + ...
            lambda_l1*sum(pressure_new) + ...
            0.5*lambda_ridge*norm(pressure_new)^2;

        % 单调重启：若加速步导致目标函数升高，则取消本次动量。
        if current_objective > objective_previous
            extrapolated_new = pressure_new;
            momentum_new = 1;
        end

        relative_change = ...
            norm(pressure_new-pressure) / ...
            max(norm(pressure_new),1);

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
        norm(G*pressure-F) / max(norm(F),eps);
    info.lipschitz_constant = lipschitz_constant;
end

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
