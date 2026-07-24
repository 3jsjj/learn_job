function number_of_nodes = ...
    nodes_per_supported_element(element_type)

    if startsWith(element_type, 'S3')
        number_of_nodes = 3;

    elseif startsWith(element_type, 'S4')
        number_of_nodes = 4;

    elseif startsWith(element_type, 'C3D4')
        number_of_nodes = 4;

    elseif startsWith(element_type, 'C3D8')
        number_of_nodes = 8;

    else
        number_of_nodes = 0;
    end
end