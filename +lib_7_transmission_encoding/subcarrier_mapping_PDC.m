%
% In 7.6.11 you are referred to 5.2.5
%
function [transmit_streams] = subcarrier_mapping_PDC(transmit_streams, physical_resource_mapping_PDC_cell, y_PDC)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    N_TS = numel(transmit_streams);
    N_b_DFT = size(cell2mat(transmit_streams(1)),1);

    % get ofdm symbol indices, last element, same for each spatial stream
    l_vec = cell2mat(physical_resource_mapping_PDC_cell(end-1));

    % loop over each spatial stream
    consumed = 0;
    for i=1:1:N_TS

        some_transmit_stream = cell2mat(transmit_streams(i));
        
        % extract PDC values for this particular transmit stream
        x_PDC_stream = cell2mat(y_PDC(i));

        % loop over each relevant ofdm symbol
        consumed_stream = 0;
        for j=1:1:numel(l_vec)

            % extract
            l = l_vec(j);
            k_i = cell2mat(physical_resource_mapping_PDC_cell(j));
            n_subc = numel(k_i);
            values = x_PDC_stream(consumed_stream+1:consumed_stream+n_subc);

            % convert
            k_i_matlab = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i);

            % map
            some_transmit_stream(k_i_matlab, l + MATLAB_INDEX_SHIFT) = values;
            consumed_stream = consumed_stream + n_subc;
        end
        
        % sanity check
        if consumed_stream ~= numel(x_PDC_stream)
            error('Not all symbols of PDC used.');
        end        
        
        % add used samples
        consumed = consumed + consumed_stream;

        % replace
        transmit_streams(i) = {some_transmit_stream};
    end

    % sanity check
    if consumed ~= numel(cell2mat(y_PDC))
        error('Not all symbols of PDC used.');
    end
end

