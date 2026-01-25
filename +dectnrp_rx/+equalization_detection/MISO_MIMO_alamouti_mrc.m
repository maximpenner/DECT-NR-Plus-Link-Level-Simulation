function [x_PDC_rev_all_rx] = MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ...
                                                     ch_estim, ...
                                                     N_RX, ...
                                                     N_eff_TX, ...
                                                     physical_resource_mapping_PDC_cell)           

    % load indices so we know which antennas were used (N_TS = N_eff_TX for tdc)
    [~, idx, prefactor] = dectnrp_6_generic_procedures.Transmit_diversity_precoding_Y(N_eff_TX);
    
    % always 98 elements in case of PCC
    n_x_PxC = numel(physical_resource_mapping_PDC_cell{end});
    
    assert(mod(n_x_PxC,2) == 0);
    
    % extracted values across all antennas
    H_all = zeros(2*N_RX, 2, n_x_PxC/2);
    y_all = zeros(2*N_RX, 1, n_x_PxC/2);
    
    for q=1:1:N_RX  
        
        % extract received samples and channel estimations for each tx antenna
        x_PDC_rev_orig = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(antenna_streams_mapped_rev(q), physical_resource_mapping_PDC_cell);
        ch_estim_PDC = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(ch_estim(q), physical_resource_mapping_PDC_cell);

        % repeat idx, see modulo operation in Table 6.3.3.2-1 to 6.3.3.2-3
        idx_rep = repmat(idx, ceil(n_x_PxC/size(idx,1)),1);
        idx_rep = idx_rep(1:n_x_PxC,:);

        % The matrix 'ch_estim_PDC' contains all channel estimates.
        % We need to extract the values at those antennas that were used.
        idx_lin = (idx_rep-1)*n_x_PxC;
        idx_lin = idx_lin + repmat((1:1:n_x_PxC)',1,2);
        h = [ch_estim_PDC(idx_lin(:,1)) ch_estim_PDC(idx_lin(:,2))];
        
        % our vector h now contains the following coefficients (notation is h_subcarrier_TXant)
        %
        %   h =     h11 h12   -> first symbol pair
        %           h21 h22   -> first symbol pair
        %           h11 h12   -> second symbol pair
        %           h21 h22   -> second symbol pair
        %           ...
        h11 = h(1:2:end,1);
        h21 = h(2:2:end,1);
        h12 = h(1:2:end,2);
        h22 = h(2:2:end,2);
        
        % Alamouti-style equations per symbol pair
        %
        % encoding
        %
        %   y1 = h11 * x1  + h12 * (-x2*)
        %   y2 = h22 * x1* + h21 * x2
        %
        %   H   = [  h11  -h12;
        %            h22*  h21*]
        %
        %   [y1; y2*] = H * [x1; x2*]
        %
        % decoding
        %
        %   H^H = [  h11*  h22;
        %           -h12*  h21]
        %
        %   [x1; x2*] = H^H * [y1; y2*]
        
        % construct H matrix for this particular antenna
        H = zeros(2,2,n_x_PxC/2);
        H(1,1,:) = reshape(h11,1,1,[]);
        H(2,1,:) = reshape(conj(h22),1,1,[]);
        H(1,2,:) = reshape(-h12,1,1,[]);
        H(2,2,:) = reshape(conj(h21),1,1,[]);
        
        % next we need y1 and y2*
        x_PDC_rev = reshape(x_PDC_rev_orig,2,[]);
        x_PDC_rev(2,:) = conj(x_PDC_rev(2,:));
        y = reshape(x_PDC_rev,2,1,[]);
        
        % fill H matrix
        H_all(q,:,:) = H(1,:,:);
        H_all(q+N_RX,:,:) = H(2,:,:);
        y_all(q,:,:) = y(1,:,:);
        y_all(q+N_RX,:,:) = y(2,:,:);
    end
    
    % output
    x_PDC_rev_all_rx = zeros(n_x_PxC,1);    
    
    % for each pair of complex symbols
    for l=1:1:(n_x_PxC/2)
        H = H_all(:,:,l);
        H_H = H';
        y = y_all(:,:,l);

        a = H_H*H;
        a = inv(a);
        x = a*H_H * y;

        x_PDC_rev_all_rx((l-1)*2+1,1) = x(1);
        x_PDC_rev_all_rx((l-1)*2+2,1) = conj(x(2));
    end

    % lastly, we need to scale with the prefactor
    x_PDC_rev_all_rx = x_PDC_rev_all_rx./prefactor;
end
