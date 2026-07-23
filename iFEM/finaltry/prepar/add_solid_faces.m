function add_solid_faces( ...
    connectivity, local_faces, ...
    node_coords, id_to_row, ...
    face_count, face_oriented)
% 将实体单元各面加入计数器。
% 两个实体共享的面计数为 2，最后会被删除。

    element_rows = zeros(1, numel(connectivity));

    for j = 1:numel(connectivity)
        element_rows(j) = id_to_row(connectivity(j));
    end

    element_centroid = mean( ...
        node_coords(element_rows, :), 1);

    for face_index = 1:numel(local_faces)
        face_ids = connectivity(local_faces{face_index});

        face_rows = zeros(1, numel(face_ids));

        for j = 1:numel(face_ids)
            face_rows(j) = id_to_row(face_ids(j));
        end

        face_points = node_coords(face_rows, :);
        face_centroid = mean(face_points, 1);

        face_normal = cross( ...
            face_points(2, :) - face_points(1, :), ...
            face_points(3, :) - face_points(1, :));

        % 将实体外表面法向调整为从单元内部指向外部
        if dot( ...
                face_normal, ...
                face_centroid - element_centroid) < 0

            face_ids = face_ids([1, end:-1:2]);
        end

        sorted_ids = sort(face_ids);
        key = sprintf('%.0f_', sorted_ids);

        if isKey(face_count, key)
            face_count(key) = face_count(key) + 1;
        else
            face_count(key) = 1;
            face_oriented(key) = face_ids;
        end
    end
end