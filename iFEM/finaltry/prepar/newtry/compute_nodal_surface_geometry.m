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