function [noise_power] = noise_estimation(antenna_streams_mapped_rev, ...
                               physical_resource_mapping_DRS_cell, ...
                               ch_estim, ...
                               N_RX, ...
                               N_eff_TX)
% we can use the channel estimate to estimate the noise power. This
% estimation is needed for the MIMO feedback as suggested in https://ieeexplore.ieee.org/document/8059870

    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    %noise
    noise_power = 0;

    % for each RX antenna
    for i=1:1:N_RX

        % create empty container
        ch_estim_i = ch_estim{i};

        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % for each transmit stream
        for j=1:1:N_eff_TX

            ch_estim_i_j = ch_estim_i(:,:,j);

            %% extract send values
            % linear indices of DRS
            linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted DRS
            values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received DRS
            y = transmit_streams_rev_i(ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab));

            %% noise estimation at pilot position
            H_hat = ch_estim_i_j(ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab));
            y_hat = H_hat.* repmat(values,1, size(y,2));

            err = y - y_hat;

            noise_power = noise_power + sum((abs(err).^2),"all")/ numel(err);
        end
    end
end
