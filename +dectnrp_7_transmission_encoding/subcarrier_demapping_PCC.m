function [x_PCC_rev] = subcarrier_demapping_PCC(any_cell_with_single_matrix, physical_resource_mapping_PCC_cell)

    % convert to matrix
    any_mat = cell2mat(any_cell_with_single_matrix);
    
    % x y coordinates
    linear_indices_matlab = cell2mat(physical_resource_mapping_PCC_cell(end));
    
    [x,y,z] = size(any_mat);

    % extract
    if z==1
        x_PCC_rev = any_mat(linear_indices_matlab);
    else
        any_mat_2d = reshape(any_mat,[x*y, z]);
        x_PCC_rev = any_mat_2d(linear_indices_matlab,:);
    end    
end
