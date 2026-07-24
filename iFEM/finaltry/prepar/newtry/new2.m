% =========================================================================
% 曲面位移插值（纯位移版本）
%
% 功能：
%   1) 从 INP 读取节点坐标和外表面三角形；
%   2) 从 CSV 读取少量测点的 NodeID、U1、U2、U3；
%   3) 将固定边界节点作为零位移锚点（可关闭）；
%   4) 沿有限元表面网格进行拉普拉斯调和插值；
%   5) 输出完整外表面的 U1、U2、U3，并绘制位移云图；
%   6) 若存在完整 Abaqus 位移文件，则自动进行逐节点误差比较。
%
% 本程序不读取 MTX，也不建立刚度矩阵。
% =========================================================================

clear;
clc;

%% 1. 用户设置

inp_filepath = 'Zhimian_addnodes_cut.inp';

% 至少包含四列：NodeID, U1, U2, U3。
% 若后面还有其他列，本程序会自动忽略。
measured_filepath = 'Abaqus_Nodal_U.csv';

% 可选：固定边界 NodeID 文件。固定节点作为 U1=U2=U3=0 的插值锚点。
fixed_filepath = 'fixed_nodes.csv';
use_fixed_zero_anchors = true;

% 可选：完整 Abaqus 表面位移，用于检查插值结果。
% 至少包含四列：NodeID, U1, U2, U3。
reference_filepath = 'Abaqus_All_Surface_Displacement.csv';
compare_with_reference_if_available = true;

% INP/HyperMesh 坐标到 ODB 全局坐标的转换。
R_inp_to_odb = diag([1, -1, -1]);

% 网格边权：相邻节点距离越短，耦合越强。
edge_weight_power = 1.0;

% 极小数值稳定项。通常无需修改。
regularization_factor = 1e-12;

% 位移显示单位。
displacement_unit = 'mm';

% 是否保存云图。
save_figures = true;

%% 2. 读取模型表面网格

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

valid_surface_triangles = ...
    all(surface_triangles_local > 0, 2);

surface_triangles_local = ...
    surface_triangles_local(valid_surface_triangles, :);

if isempty(surface_triangles_local)
    error('未找到可用的外表面三角形。');
end

fprintf('读取模型节点：%d\n', numel(node_ids));
fprintf('外表面节点：%d\n', N_surface_nodes);
fprintf('外表面三角形：%d\n', size(surface_triangles_local, 1));

%% 3. 读取测点位移

measured_data = readmatrix(measured_filepath);

if size(measured_data, 2) < 4
    error('测点文件至少需要四列：NodeID, U1, U2, U3。');
end

measured_ids_raw = measured_data(:, 1);
U_measured_raw = measured_data(:, 2:4);

valid_measured_rows = ...
    isfinite(measured_ids_raw) & ...
    all(isfinite(U_measured_raw), 2);

measured_ids_raw = round(measured_ids_raw(valid_measured_rows));
U_measured_raw = U_measured_raw(valid_measured_rows, :);

% 同一 NodeID 若重复出现，则对位移取平均值。
[measured_ids, ~, measured_group] = unique(measured_ids_raw);
U_measured = zeros(numel(measured_ids), 3);

for component_index = 1:3
    U_measured(:, component_index) = accumarray( ...
        measured_group, ...
        U_measured_raw(:, component_index), ...
        [numel(measured_ids), 1], ...
        @mean);
end

if numel(measured_ids) < numel(measured_ids_raw)
    warning('测点文件存在重复 NodeID，已对重复位移取平均值。');
end

% NodeID -> INP节点行号
[measured_exists_in_model, measured_node_rows] = ...
    ismember(measured_ids, node_ids);

if any(~measured_exists_in_model)
    warning('%d 个测点 NodeID 不在 INP 中，已剔除。', ...
        nnz(~measured_exists_in_model));
end

measured_ids = measured_ids(measured_exists_in_model);
measured_node_rows = measured_node_rows(measured_exists_in_model);
U_measured = U_measured(measured_exists_in_model, :);

% 只保留外表面测点
[measured_is_on_surface, measured_surface_rows] = ...
    ismember(measured_node_rows, surface_node_rows);

if any(~measured_is_on_surface)
    warning('%d 个测点不属于识别出的外表面，已剔除。', ...
        nnz(~measured_is_on_surface));
end

measured_ids = measured_ids(measured_is_on_surface);
measured_surface_rows = measured_surface_rows(measured_is_on_surface);
U_measured = U_measured(measured_is_on_surface, :);

if isempty(measured_surface_rows)
    error('没有可用于插值的外表面测点。');
end

fprintf('有效位移测点：%d\n', numel(measured_surface_rows));

%% 4. 建立已知位移锚点

known_surface_rows = measured_surface_rows(:);
known_displacements = U_measured;
known_source = ones(numel(known_surface_rows), 1);  % 1=测量点，2=固定点

if use_fixed_zero_anchors && isfile(fixed_filepath)

    fixed_data = readmatrix(fixed_filepath);
    fixed_ids = unique(round(fixed_data(isfinite(fixed_data))));
    fixed_ids = fixed_ids(fixed_ids > 0);

    [fixed_exists_in_model, fixed_node_rows] = ...
        ismember(fixed_ids, node_ids);

    if any(~fixed_exists_in_model)
        warning('%d 个固定节点 NodeID 不在 INP 中，已剔除。', ...
            nnz(~fixed_exists_in_model));
    end

    fixed_ids = fixed_ids(fixed_exists_in_model);
    fixed_node_rows = fixed_node_rows(fixed_exists_in_model);

    [fixed_is_on_surface, fixed_surface_rows] = ...
        ismember(fixed_node_rows, surface_node_rows);

    fixed_ids = fixed_ids(fixed_is_on_surface);
    fixed_surface_rows = fixed_surface_rows(fixed_is_on_surface);

    % 测量值优先；只添加尚未作为测点的固定节点。
    fixed_not_measured = ...
        ~ismember(fixed_surface_rows, known_surface_rows);

    fixed_surface_rows_to_add = ...
        fixed_surface_rows(fixed_not_measured);

    known_surface_rows = [ ...
        known_surface_rows; ...
        fixed_surface_rows_to_add(:)];

    known_displacements = [ ...
        known_displacements; ...
        zeros(numel(fixed_surface_rows_to_add), 3)];

    known_source = [ ...
        known_source; ...
        2 * ones(numel(fixed_surface_rows_to_add), 1)];

    % 如果固定节点同时出现在测点文件中，但测量值明显非零，给出提示。
    [fixed_also_measured, measured_location] = ...
        ismember(fixed_surface_rows, measured_surface_rows);

    if any(fixed_also_measured)
        fixed_measured_u = ...
            U_measured(measured_location(fixed_also_measured), :);

        nonzero_fixed_measurements = ...
            vecnorm(fixed_measured_u, 2, 2) > 1e-10;

        if any(nonzero_fixed_measurements)
            warning(['%d 个固定节点在测点文件中的位移不为零；', ...
                '当前保留测量值，不强制改为零。'], ...
                nnz(nonzero_fixed_measurements));
        end
    end

    fprintf('新增固定零位移锚点：%d\n', ...
        numel(fixed_surface_rows_to_add));

elseif use_fixed_zero_anchors
    warning('未找到 fixed_nodes.csv，当前只使用测点进行插值。');
end

%% 5. 构造表面网格拉普拉斯矩阵

tri = surface_triangles_local;

edge_list = [ ...
    tri(:, [1, 2]); ...
    tri(:, [2, 3]); ...
    tri(:, [3, 1])];

edge_list = sort(edge_list, 2);
edge_list = unique(edge_list, 'rows');

edge_vector = ...
    surface_coords(edge_list(:,1), :) - ...
    surface_coords(edge_list(:,2), :);

edge_length = vecnorm(edge_vector, 2, 2);

positive_edge_length = edge_length(edge_length > 0);

if isempty(positive_edge_length)
    error('表面网格边长度无效。');
end

minimum_edge_length = max( ...
    min(positive_edge_length) * 1e-12, ...
    eps);

edge_weight = 1 ./ ...
    max(edge_length, minimum_edge_length).^edge_weight_power;

W = sparse( ...
    [edge_list(:,1); edge_list(:,2)], ...
    [edge_list(:,2); edge_list(:,1)], ...
    [edge_weight; edge_weight], ...
    N_surface_nodes, ...
    N_surface_nodes);

L = spdiags(sum(W, 2), 0, ...
    N_surface_nodes, N_surface_nodes) - W;

%% 6. 调和插值 U1、U2、U3

U_surface_3D = nan(N_surface_nodes, 3);
U_surface_3D(known_surface_rows, :) = known_displacements;

is_known = false(N_surface_nodes, 1);
is_known(known_surface_rows) = true;

% 分连通区域处理，避免孤立表面导致奇异矩阵。
mesh_graph = graph(W);
component_id = conncomp(mesh_graph)';
component_count = max(component_id);

fallback_node_count = 0;

for current_component = 1:component_count

    component_nodes = find(component_id == current_component);
    component_known = component_nodes(is_known(component_nodes));
    component_unknown = component_nodes(~is_known(component_nodes));

    if isempty(component_unknown)
        continue;
    end

    if isempty(component_known)
        % 当前连通区域没有锚点：使用空间最近的全局锚点值。
        for local_index = 1:numel(component_unknown)
            current_node = component_unknown(local_index);

            coordinate_difference = bsxfun( ...
                @minus, ...
                surface_coords(known_surface_rows, :), ...
                surface_coords(current_node, :));

            squared_distance = sum(coordinate_difference.^2, 2);
            [~, nearest_known_index] = min(squared_distance);

            U_surface_3D(current_node, :) = ...
                known_displacements(nearest_known_index, :);
        end

        fallback_node_count = ...
            fallback_node_count + numel(component_unknown);

        continue;
    end

    A = L(component_unknown, component_unknown);
    B = -L(component_unknown, component_known) * ...
        U_surface_3D(component_known, :);

    diagonal_mean = full(mean(diag(A)));

    if ~isfinite(diagonal_mean) || diagonal_mean <= 0
        diagonal_mean = 1;
    end

    A = A + ...
        regularization_factor * diagonal_mean * ...
        speye(size(A, 1));

    U_surface_3D(component_unknown, :) = A \ B;
end

if fallback_node_count > 0
    warning(['%d 个节点所在的连通表面没有任何位移锚点，', ...
        '已使用最近锚点值。'], fallback_node_count);
end

if any(~isfinite(U_surface_3D(:)))
    error('位移插值后仍存在无效值，请检查表面网格连通性。');
end

% 确保测点值完全保持不变。
U_surface_3D(measured_surface_rows, :) = U_measured;

fprintf('完成插值节点：%d\n', ...
    N_surface_nodes - numel(measured_surface_rows));

%% 7. 输出完整外表面位移

is_measured_node = false(N_surface_nodes, 1);
is_measured_node(measured_surface_rows) = true;

is_fixed_anchor = false(N_surface_nodes, 1);
is_fixed_anchor(known_surface_rows(known_source == 2)) = true;

U_magnitude = vecnorm(U_surface_3D, 2, 2);

result_table = table( ...
    surface_node_ids(:), ...
    surface_coords(:,1), ...
    surface_coords(:,2), ...
    surface_coords(:,3), ...
    U_surface_3D(:,1), ...
    U_surface_3D(:,2), ...
    U_surface_3D(:,3), ...
    U_magnitude, ...
    is_measured_node, ...
    is_fixed_anchor, ...
    'VariableNames', { ...
        'NodeID', ...
        'X', 'Y', 'Z', ...
        'U1_Interpolated', ...
        'U2_Interpolated', ...
        'U3_Interpolated', ...
        'U_Magnitude', ...
        'IsMeasuredNode', ...
        'IsFixedAnchor'});

output_filepath = 'Surface_Displacement_Interpolation.csv';
writetable(result_table, output_filepath);

fprintf('结果已保存：%s\n', output_filepath);

%% 8. 绘制 U1、U2、U3 和位移模长云图

plot_values = { ...
    U_surface_3D(:,1), ...
    U_surface_3D(:,2), ...
    U_surface_3D(:,3), ...
    U_magnitude};

plot_names = {'U1', 'U2', 'U3', 'UMagnitude'};
plot_titles = { ...
    '表面插值位移 U1', ...
    '表面插值位移 U2', ...
    '表面插值位移 U3', ...
    '表面插值位移模长'};

for plot_index = 1:numel(plot_values)

    current_value = plot_values{plot_index};

    figure( ...
        'Name', plot_titles{plot_index}, ...
        'Color', 'w', ...
        'Position', [100, 100, 900, 700]);

    trisurf( ...
        surface_triangles_local, ...
        surface_coords(:,1), ...
        surface_coords(:,2), ...
        surface_coords(:,3), ...
        current_value, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none');

    hold on;

    scatter3( ...
        surface_coords(measured_surface_rows,1), ...
        surface_coords(measured_surface_rows,2), ...
        surface_coords(measured_surface_rows,3), ...
        30, 'k', 'filled');

    axis equal;
    axis tight;
    grid on;
    box on;
    view(45, 30);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(plot_titles{plot_index});

    cb = colorbar;
    ylabel(cb, sprintf('%s (%s)', ...
        plot_names{plot_index}, displacement_unit));

    colormap(jet(256));

    if plot_index <= 3
        current_limit = max(abs(current_value));
        if current_limit > 0
            caxis([-current_limit, current_limit]);
        end
    end

    legend('插值曲面', '测量节点', ...
        'Location', 'best');

    hold off;

    if save_figures
        output_image_name = sprintf( ...
            'Interpolated_%s_Cloud.png', ...
            plot_names{plot_index});

        exportgraphics(gcf, output_image_name, ...
            'Resolution', 300);

        fprintf('已保存：%s\n', output_image_name);
    end
end

%% 9. 可选：与完整 Abaqus 表面位移比较

if compare_with_reference_if_available && ...
        isfile(reference_filepath)

    reference_data = readmatrix(reference_filepath);

    if size(reference_data, 2) < 4
        warning('完整位移参考文件少于四列，已跳过比较。');
    else
        reference_ids_raw = reference_data(:,1);
        U_reference_raw = reference_data(:,2:4);

        valid_reference_rows = ...
            isfinite(reference_ids_raw) & ...
            all(isfinite(U_reference_raw), 2);

        reference_ids = ...
            round(reference_ids_raw(valid_reference_rows));
        U_reference = ...
            U_reference_raw(valid_reference_rows, :);

        [matched_surface, reference_location] = ...
            ismember(surface_node_ids, reference_ids);

        compare_rows = find(matched_surface);
        U_interpolated_compare = U_surface_3D(compare_rows, :);
        U_reference_compare = ...
            U_reference(reference_location(compare_rows), :);

        interpolation_error = ...
            U_interpolated_compare - U_reference_compare;

        fprintf('\n');
        fprintf('========== 完整位移场比较 ==========\n');
        fprintf('匹配节点数量：%d / %d\n', ...
            numel(compare_rows), N_surface_nodes);

        component_names = {'U1', 'U2', 'U3'};

        for component_index = 1:3
            relative_error = ...
                norm(interpolation_error(:,component_index)) / ...
                max(norm(U_reference_compare(:,component_index)), eps);

            current_correlation = corrcoef( ...
                U_interpolated_compare(:,component_index), ...
                U_reference_compare(:,component_index));

            if size(current_correlation, 1) == 2
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
            norm(interpolation_error, 'fro') / ...
            max(norm(U_reference_compare, 'fro'), eps);

        overall_cosine_similarity = ...
            dot(U_interpolated_compare(:), ...
                U_reference_compare(:)) / ...
            max(norm(U_interpolated_compare(:)) * ...
                norm(U_reference_compare(:)), eps);

        fprintf('整体相对L2误差：%.2f %%\n', ...
            100 * overall_relative_error);
        fprintf('整体余弦相似度：%.6f\n', ...
            overall_cosine_similarity);
        fprintf('====================================\n');

        reference_error_magnitude = ...
            vecnorm(interpolation_error, 2, 2);

        comparison_table = table( ...
            surface_node_ids(compare_rows), ...
            U_reference_compare(:,1), ...
            U_interpolated_compare(:,1), ...
            interpolation_error(:,1), ...
            U_reference_compare(:,2), ...
            U_interpolated_compare(:,2), ...
            interpolation_error(:,2), ...
            U_reference_compare(:,3), ...
            U_interpolated_compare(:,3), ...
            interpolation_error(:,3), ...
            reference_error_magnitude, ...
            is_measured_node(compare_rows), ...
            'VariableNames', { ...
                'NodeID', ...
                'U1_Abaqus', 'U1_Interpolated', 'U1_Error', ...
                'U2_Abaqus', 'U2_Interpolated', 'U2_Error', ...
                'U3_Abaqus', 'U3_Interpolated', 'U3_Error', ...
                'ErrorMagnitude', ...
                'IsMeasuredNode'});

        writetable(comparison_table, ...
            'Abaqus_vs_Interpolated_Displacement.csv');

        fprintf('比较结果已保存：');
        fprintf('Abaqus_vs_Interpolated_Displacement.csv\n');
    end
end

fprintf('\n程序运行完成。\n');
