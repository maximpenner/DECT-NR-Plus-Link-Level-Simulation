%
% In 7.5.9 you are referred to 5.2.4
%
function [transmit_streams] = subcarrier_mapping_PCC(transmit_streams, physical_resource_mapping_PCC_cell, y_PCC)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    N_TS = numel(transmit_streams);
    N_b_DFT = size(cell2mat(transmit_streams(1)),1);    

    % get all relevant ofdm symbols, equal for each transmit stream
    l_vec = cell2mat(physical_resource_mapping_PCC_cell(end-1));
    
    % map for each transmit stream
    for i=1:1:N_TS
        
        some_transmit_stream = cell2mat(transmit_streams(i));
        
        % extract values for this transmit stream
        y_PCC_this_stream = cell2mat(y_PCC(i));        

        % loop over each relevant ofdm symbol
        for j=1:1:numel(l_vec)

            % extract
            l = l_vec(j);
            k_i = cell2mat(physical_resource_mapping_PCC_cell(j));
            n_subc = numel(k_i);
            values = y_PCC_this_stream(1:n_subc);
            y_PCC_this_stream(1:n_subc) = [];

            % convert
            k_i_matlab = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i);

            % map
            some_transmit_stream(k_i_matlab, l + MATLAB_INDEX_SHIFT) = values;
        end
        
        % replace
        transmit_streams(i) = {some_transmit_stream};        

        % sanity check
        if numel(y_PCC_this_stream) ~= 0
            error('Not all symbols of PCC used.');
        end
    end
end

