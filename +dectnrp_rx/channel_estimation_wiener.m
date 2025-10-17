function [ch_estim] = channel_estimation_wiener(antenna_streams_mapped_rev, ...
                                                physical_resource_mapping_DRS_cell, ...
                                                wiener, ...
                                                N_RX, ...
                                                N_eff_TX)

    % we need the size of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    % output
    ch_estim = cell(N_RX,1);

    % for each rx antenna
    for i=1:1:N_RX

        % create empty container
        ch_estim_i = zeros(N_b_DFT, N_PACKET_symb, N_eff_TX);

        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % for each tx antenna
        for j=1:1:N_eff_TX
            %% extract send values
            % linear indices of drs
            linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted drs
            values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received drs
            y = transmit_streams_rev_i(ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab));

            %% least squares estimation at the pilot positions
            ls = y./repmat(values,1, size(y,2));
            ls = reshape(ls,[],1);
            
            %% create on big matrix of the estimated channel at the pilot positions
            ls = reshape(ls,1,1,numel(ls));
            ls = repmat(ls,N_b_DFT, N_PACKET_symb, 1);

            %% wiener interpolate
            ch_estim_i(:,:,j) = sum(ls.*cell2mat(wiener(j)),3);
        end
        ch_estim(i) = {ch_estim_i};
    end
end
