function [node_ids, node_coords] = get_model(inp_filepath)
    % 从 Abaqus .inp 文件中高效提取节点坐标
    fid = fopen(inp_filepath, 'r');
    if fid == -1
        error('无法打开指定的 .inp 文件');
    end

    node_data = [];
    is_node_section = false;

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        % 寻找 *Node 关键字
        if startsWith(line, '*Node', 'IgnoreCase', true)
            is_node_section = true;
            continue;
        elseif startsWith(line, '*') && is_node_section
            % 如果遇到其他关键字（以 * 开头），说明节点数据块结束了
            is_node_section = false;
        end

        if is_node_section && ~isempty(line)
            % 解析当前行的数字（格式为：ID, X, Y, Z）
            vals = str2num(line); %#ok<ST2NM>
            if ~isempty(vals)
                node_data = [node_data; vals]; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    node_ids = node_data(:, 1);
    node_coords = node_data(:, 2:end); % 得到 [N_nodes, 3] 矩阵
end