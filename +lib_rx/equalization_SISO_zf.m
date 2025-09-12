function [x_PDC_rev] = equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell)

    % sanity check
    if numel(antenna_streams_mapped_rev) ~= 1 || numel(ch_estim) ~= 1
        error('Must be SISO.');
    end
    
    % extract
    antenna_streams_mapped_rev_mat = cell2mat(antenna_streams_mapped_rev);
    ch_estim_mat = cell2mat(ch_estim);

    % zf
    transmit_streams_rev = {antenna_streams_mapped_rev_mat./ch_estim_mat};
    
    % extract PDC
    x_PDC_rev = lib_7_transmission_encoding.subcarrier_demapping_PDC(transmit_streams_rev, physical_resource_mapping_PDC_cell);
end

