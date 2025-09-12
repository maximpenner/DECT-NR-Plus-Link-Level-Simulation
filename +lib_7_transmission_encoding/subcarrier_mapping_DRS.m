%
% 7.4 and 5.2.3
%
function [transmit_streams] = subcarrier_mapping_DRS(transmit_streams, physical_resource_mapping_DRS_cell)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    N_TS = numel(transmit_streams);
    N_b_DFT = size(cell2mat(transmit_streams(1)),1);

    % mapping different for each transmit stream
    for i=1:1:N_TS

        some_transmit_stream = cell2mat(transmit_streams(i));

        k_i_mat = cell2mat(physical_resource_mapping_DRS_cell(i,1));
        l_vec = cell2mat(physical_resource_mapping_DRS_cell(i,2));

        % stays the same in a transmit stream
        values = cell2mat(physical_resource_mapping_DRS_cell(i,3));

        for j=1:1:numel(l_vec)
            
            % extract
            l = l_vec(j);
            k_i = k_i_mat(:,j);

            % convert
            k_i_column_vec_matlab = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i);

            % map
            some_transmit_stream(k_i_column_vec_matlab, l + MATLAB_INDEX_SHIFT) = values;
        end

        % replace
        transmit_streams(i) = {some_transmit_stream};
    end
end

