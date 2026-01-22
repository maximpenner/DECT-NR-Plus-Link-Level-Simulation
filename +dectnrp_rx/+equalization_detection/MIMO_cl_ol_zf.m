function [x_PDC_rev] = MIMO_cl_ol_zf(antenna_streams_mapped_rev, ...
                                                     ch_estim, ...
                                                     N_RX, ...
                                                     N_eff_TX, ...
                                                     physical_resource_mapping_PDC_cell, ...
                                                     N_ss)           

    assert(N_RX >= N_eff_TX, "number of transmit and receive antennas must be equal for open loop MIMO!");

    % always 98 elements in case of PCC
    n_x_PxC = numel(physical_resource_mapping_PDC_cell{end});
    
    assert(mod(n_x_PxC,2) == 0);
    
    % extracted values across all antennas
    H_all = zeros(N_RX, N_eff_TX ,n_x_PxC);
    y_all = zeros(N_RX, n_x_PxC);
    
    for q=1:1:N_RX  
        
        % extract received samples and channel estimations for each tx antenna
        x_PDC_rev_orig = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(antenna_streams_mapped_rev(q), physical_resource_mapping_PDC_cell);
        ch_estim_PDC = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(ch_estim(q), physical_resource_mapping_PDC_cell);

        % store received symbols
        y_all(q, :) = x_PDC_rev_orig(:).';

        % extract channel
        for t = 1:N_eff_TX
            H_all(q,t,:) = ch_estim_PDC(:,t);
        end
    end
    
    % output
    x_PDC_rev_all_rx = zeros(N_eff_TX,n_x_PxC);    
    
    % for each subcarrier
    for k=1:1:(n_x_PxC)
        H_k = squeeze(H_all(:,:,k));
        y_k = y_all(:,k);

        % zero forcing (solve linear equation or inv)
        % x_PDC_rev_all_rx(:, k) = (H_k' * H_k) \ (H_k' * y_k);
        H_tH_inv = inv(H_k' * H_k);
        x_PDC_rev_all_rx(:, k) = H_tH_inv * (H_k' * y_k);
    end

    % Interleave according to
    % +dectnrp_6_generic_procedures/Spatial_Multiplexing

    x_PDC_rev = zeros(n_x_PxC * N_ss, 1);
    for i = 1:N_ss
        x_PDC_rev(i:N_ss:end) = x_PDC_rev_all_rx(i, :);
    end
        
end
