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