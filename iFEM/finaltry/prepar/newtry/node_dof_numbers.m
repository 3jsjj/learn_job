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
