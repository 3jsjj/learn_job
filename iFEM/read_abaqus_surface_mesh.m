function [node_ids, node_coords, surface_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath)
% 读取常见 Abaqus 节点和单元，并提取表面三角形。
%
% 当前支持：
%   壳单元：S3、S3R、S4、S4R
%   实体单元：C3D4、C3D8、C3D8R
%
% 对实体单元：
%   自动删除两个实体共享的内部面，只保留外表面。
%
% surface_triangles 中保存的是 node_coords 的行号，而不是 NodeID。

    raw_text = fileread(inp_filepath);
    lines = splitlines(raw_text);

    node_ids = zeros(0, 1);
    node_coords = zeros(0, 3);

    element_types = cell(0, 1);
    element_connectivity = cell(0, 1);

    current_mode = '';
    current_element_type = '';
    expected_element_nodes = 0;
    element_buffer = [];

    for line_index = 1:numel(lines)
        line = strtrim(lines{line_index});

        if isempty(line) || startsWith(line, '**')
            continue;
        end

        if startsWith(line, '*')
            upper_line = upper(line);

            if startsWith(upper_line, '*NODE')
                current_mode = 'node';
                current_element_type = '';
                expected_element_nodes = 0;
                element_buffer = [];

            elseif startsWith(upper_line, '*ELEMENT')
                current_mode = 'element';

                token = regexp( ...
                    upper_line, ...
                    'TYPE\s*=\s*([^,\s]+)', ...
                    'tokens', 'once');

                if isempty(token)
                    current_element_type = '';
                    expected_element_nodes = 0;
                else
                    current_element_type = token{1};
                    expected_element_nodes = ...
                        nodes_per_supported_element( ...
                            current_element_type);
                end

                element_buffer = [];
            else
                current_mode = '';
                current_element_type = '';
                expected_element_nodes = 0;
                element_buffer = [];
            end

            continue;
        end

        values = sscanf( ...
            strrep(line, ',', ' '), '%f')';

        if isempty(values)
            continue;
        end

        if strcmp(current_mode, 'node')
            if numel(values) < 3
                continue;
            end

            node_ids(end + 1, 1) = values(1);

            xyz = zeros(1, 3);
            number_of_coordinates = min(numel(values) - 1, 3);
            xyz(1:number_of_coordinates) = ...
                values(2:1 + number_of_coordinates);

            node_coords(end + 1, :) = xyz;

        elseif strcmp(current_mode, 'element') && ...
                expected_element_nodes > 0

            element_buffer = [element_buffer, values]; %#ok<AGROW>

            values_per_element = ...
                expected_element_nodes + 1;

            while numel(element_buffer) >= values_per_element
                one_element = ...
                    element_buffer(1:values_per_element);

                element_buffer(1:values_per_element) = [];

                element_types{end + 1, 1} = ...
                    current_element_type;

                element_connectivity{end + 1, 1} = ...
                    one_element(2:end);
            end
        end
    end

    if isempty(node_ids)
        error('未能从 inp 文件读取节点。');
    end

    if isempty(element_connectivity)
        error(['未读取到受支持的单元。', newline, ...
               '当前支持 S3/S4/C3D4/C3D8。']);
    end

    id_to_row = containers.Map( ...
        'KeyType', 'double', ...
        'ValueType', 'double');

    for n = 1:numel(node_ids)
        id_to_row(node_ids(n)) = n;
    end

    shell_faces = cell(0, 1);

    face_count = containers.Map( ...
        'KeyType', 'char', ...
        'ValueType', 'double');

    face_oriented = containers.Map( ...
        'KeyType', 'char', ...
        'ValueType', 'any');

    for element_index = 1:numel(element_connectivity)
        element_type = element_types{element_index};
        connectivity = element_connectivity{element_index};

        if startsWith(element_type, 'S3')
            shell_faces{end + 1, 1} = connectivity(1:3);

        elseif startsWith(element_type, 'S4')
            shell_faces{end + 1, 1} = connectivity(1:4);

        elseif startsWith(element_type, 'C3D4')
            local_faces = { ...
                [1, 2, 3], ...
                [1, 4, 2], ...
                [2, 4, 3], ...
                [3, 4, 1]};

            add_solid_faces( ...
                connectivity, local_faces, ...
                node_coords, id_to_row, ...
                face_count, face_oriented);

        elseif startsWith(element_type, 'C3D8')
            local_faces = { ...
                [1, 2, 3, 4], ...
                [5, 8, 7, 6], ...
                [1, 5, 6, 2], ...
                [2, 6, 7, 3], ...
                [3, 7, 8, 4], ...
                [4, 8, 5, 1]};

            add_solid_faces( ...
                connectivity, local_faces, ...
                node_coords, id_to_row, ...
                face_count, face_oriented);
        end
    end

    surface_faces = shell_faces;

    solid_face_keys = keys(face_count);

    for key_index = 1:numel(solid_face_keys)
        key = solid_face_keys{key_index};

        if face_count(key) == 1
            surface_faces{end + 1, 1} = ...
                face_oriented(key);
        end
    end

    surface_triangles = zeros(0, 3);

    for face_index = 1:numel(surface_faces)
        face_node_ids = surface_faces{face_index};

        face_rows = zeros(1, numel(face_node_ids));

        for j = 1:numel(face_node_ids)
            if ~isKey(id_to_row, face_node_ids(j))
                error('单元引用了不存在的节点 ID。');
            end

            face_rows(j) = id_to_row(face_node_ids(j));
        end

        if numel(face_rows) == 3
            surface_triangles(end + 1, :) = ...
                face_rows;

        elseif numel(face_rows) == 4
            surface_triangles(end + 1, :) = ...
                face_rows([1, 2, 3]);

            surface_triangles(end + 1, :) = ...
                face_rows([1, 3, 4]);
        end
    end

    if isempty(surface_triangles)
        error('没有生成有效的表面三角形。');
    end
end