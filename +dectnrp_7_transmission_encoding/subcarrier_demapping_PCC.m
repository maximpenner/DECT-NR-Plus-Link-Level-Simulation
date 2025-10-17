function [x_PCC_rev] = subcarrier_demapping_PCC(any_cell_with_single_matrix, physical_resource_mapping_PCC_cell)

    % % Technical Specification assumes first index is 0, matlab 1
    % MATLAB_INDEX_SHIFT = 1;
    % 
    % x_PCC_rev = [];
    % 
    % any_mat = cell2mat(any_cell_with_single_matrix);
    % N_b_DFT = size(any_mat,1);
    % 
    % % get OFDM symbol indices, last element
    % l_vec = cell2mat(physical_resource_mapping_PCC_cell(end-1));
    % 
    % % loop over each relevant OFDM symbol
    % for i=1:1:numel(l_vec)
    % 
    %     % extract
    %     l = l_vec(i);
    %     k_i = cell2mat(physical_resource_mapping_PCC_cell(i));
    % 
    %     % convert
    %     k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);
    % 
    %     % map
    %     x_PCC_rev = [x_PCC_rev; any_mat(k_i_matlab, l + MATLAB_INDEX_SHIFT)];
    % end

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
