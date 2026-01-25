function [x_PDC_rev] = SIMO_mrc(antenna_streams_mapped_rev, ...
                                ch_estim, ...
                                physical_resource_mapping_PDC_cell, ...
                                N_RX)

    % source: https://en.wikipedia.org/wiki/Maximal-ratio_combining
    
    % equalized and weighted received signals
    numerator_equalized_symbol = zeros(size(cell2mat(antenna_streams_mapped_rev(1))));
    
    % total power across all signals
    denominator_power = zeros(size(numerator_equalized_symbol));
  
    for i=1:1:N_RX
        ch_estim_i = cell2mat(ch_estim(i));
        antenna_streams_mapped_rev_i = cell2mat(antenna_streams_mapped_rev(i));
        
        numerator_equalized_symbol = numerator_equalized_symbol + conj(ch_estim_i).*antenna_streams_mapped_rev_i;
        denominator_power = denominator_power + abs(ch_estim_i).^2;
    end
    
    % scale with power of all channel estimates
    antenna_streams_mapped_combined_rev = {numerator_equalized_symbol./denominator_power};
    
    % extract PDC
    x_PDC_rev = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(antenna_streams_mapped_combined_rev, physical_resource_mapping_PDC_cell);    
end
