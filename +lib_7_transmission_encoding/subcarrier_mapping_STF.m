%
% 7.3 and 5.2.2
%
function [transmit_streams] = subcarrier_mapping_STF(transmit_streams, physical_resource_mapping_STF_cell)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    % STF is mapped into first transmit stream
    first_transmit_stream = cell2mat(transmit_streams(1));

    N_b_DFT = size(first_transmit_stream,1);

    % extract
    k_i = cell2mat(physical_resource_mapping_STF_cell(1));
    l = cell2mat(physical_resource_mapping_STF_cell(2));
    values = cell2mat(physical_resource_mapping_STF_cell(3));

    % convert
    k_i_matlab = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i);

    % map
    first_transmit_stream(k_i_matlab, l + MATLAB_INDEX_SHIFT) = values;

    % replace
    transmit_streams(1) = {first_transmit_stream};
end

