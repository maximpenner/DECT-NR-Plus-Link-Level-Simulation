function antenna_streams_mapped_rev_cfo = cfo_residual(antenna_streams_mapped_rev, ...
                                                       physical_resource_mapping_DRS_cell, ...
                                                       physical_resource_mapping_STF_cell, ...
                                                       N_RX, ...
                                                       N_eff_TX)
   
    % we need the dimensions of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    %% STF preparation

    % For each RX antenna, the same STF is received.
    % Get the STF indices and the sent STF complex values.
    STF_linear_indices_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, physical_resource_mapping_STF_cell{1});
    STF_values                = physical_resource_mapping_STF_cell{3};

    %% DRS preparation

    % How many OFDM symbols contain DRS symbols?
    n_DRS_symbols = numel(physical_resource_mapping_DRS_cell{1,2});

    % extract the OFDM symbols indices which contains DRS values
    DRS_symbol_indices = zeros(N_eff_TX, n_DRS_symbols);
    for i=1:1:N_eff_TX
        DRS_symbol_indices(i,:) = physical_resource_mapping_DRS_cell{i,2};
    end

    % for each TX antenna, determine the symbol to symbol distance between OFDM symbols with DRS values
    if N_eff_TX <= 2
        N_step = 5;
    else
        N_step = 10;
    end

    STF_and_DRS_symbol_2_symbol_spacing = N_step * ones(1, n_DRS_symbols);

    % Since we also include the STF for the residual CFO calculation, we need to include the distance between the STF and the first DRS symbols.
    % This distance is one OFDM symbol for N_eff_TX = 1, 2 and 4.
    % For N_eff_TX = 8, however, the average distance is 1.5 OFDM symbols (one symbol for transmit streams 1 to 4 and two symbols for transmit streams 5 to 8).
    if N_eff_TX < 8
        STF_and_DRS_symbol_2_symbol_spacing(1) = 1;
    else
        STF_and_DRS_symbol_2_symbol_spacing(1) = 1.5;
    end

    %% output preparation

    % init CFO corrected output container
    antenna_streams_mapped_rev_cfo = cell(N_RX,1);

    %% the residual CFO is corrected for each RX antenna individually
    for i=1:1:N_RX

        % extract received frequency domain samples for RX antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        %% extract received complex STF values (sent values were extracted above)

        % received
        y_STF = transmit_streams_rev_i(STF_linear_indices_matlab);

        % equation (8) from https://openOFDM.readthedocs.io/en/latest/_downloads/vtc04_freq_offset.pdf
        weighted_STF = STF_values.*conj(y_STF);

        % this line can be used to debug the angle
        %angles_STF = get_wrapped_angles_in_deg(weighted_STF);

        % equation (8) from https://openOFDM.readthedocs.io/en/latest/_downloads/vtc04_freq_offset.pdf
        %STF_pilot_value_product = sum(weighted_STF);

        % at each rx antenna we receive DRS symbols from each tx antenna
        symbol_2_symbol_rotation_mat = zeros(N_eff_TX, n_DRS_symbols);
        for j=1:1:N_eff_TX

            %% extract sent and received DRS values

            % linear indices of drs values
            DRS_linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted drs
            DRS_values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received drs
            y_DRS = transmit_streams_rev_i(DRS_linear_indices_matlab);

            %% get the angle per symbol

            % equation (8) from https://openOFDM.readthedocs.io/en/latest/_downloads/vtc04_freq_offset.pdf
            weighted_DRS = repmat(DRS_values,1, n_DRS_symbols).*conj(y_DRS);

            % this line can be used to debug the angle
            %angles_DRS = get_wrapped_angles_in_deg(weighted_DRS);

            % equation (8) from https://openOFDM.readthedocs.io/en/latest/_downloads/vtc04_freq_offset.pdf
            %DRS_pilot_value_product = sum(weighted_DRS, 1);

            %% do no average across subcarriers just yet, instead determine the angle change in time domain
            symbol_2_symbol_rotations = [weighted_STF, weighted_DRS];
            symbol_2_symbol_rotations = symbol_2_symbol_rotations(:, 2:end).*conj(symbol_2_symbol_rotations(:, 1:end-1));

            % this line can be used to debug the angle
            %angles_s2s = get_wrapped_angles_in_deg(symbol_2_symbol_rotation);

            %% now average across all subcarriers of this spatial stream
            symbol_2_symbol_rotation_mat(j,:) = sum(symbol_2_symbol_rotations, 1);

            % this line can be used to debug the angle
            %angles_s2s = get_wrapped_angles_in_deg(symbol_2_symbol_rotation_mat);
        end

        %% average across all transmit streams
        symbol_2_symbol_rotations = sum(symbol_2_symbol_rotation_mat, 1);

        % this line can be used to debugthe angle
        %angles_s2s = get_wrapped_angles_in_deg(symbol_2_symbol_rotation);
        %angles_s2s = get_wrapped_angles_in_deg(symbol_2_symbol_rotation)./STF_and_DRS_symbol_2_symbol_spacing;

        %% we now have the CFO estimate across all transmit streams for this particular RX antenna

        % So far, we have only muluiplied and summed up complex numbers.
        % Their absolute becomes an implicit weights for the next step: determining the angle.
        symbol_2_symbol_angles = angle(symbol_2_symbol_rotations);
    
        % this is very important for random phases
        symbol_2_symbol_angles = dectnrp_util.wrap2_plus_pi_minus_pi(symbol_2_symbol_angles);

        % this line can be used to debug the angle
        %angles_s2s = rad2deg(symbol_2_symbol_angles);

        % normalize the angle to the number of OFDM symbols spacing
        symbol_2_symbol_angles_normalized = symbol_2_symbol_angles./STF_and_DRS_symbol_2_symbol_spacing;

        % this line can be used to debug the angle
        %angles_s2s = rad2deg(symbol_2_symbol_angles_normalized);
    
        %% based on the angles between STF and DRS symbols, we now calculate a derotation matrix for all OFDM symbols

        % It can happen that after the last DRS symbols there are still OFDM symbols following.
        % In this case, we have extrapolate derotation values.
        n_OFDM_symbols_after_last_DRS_symbol = N_PACKET_symb - DRS_symbol_indices(end,end) - 1;

        if N_eff_TX == 8
            n_OFDM_symbols_after_last_DRS_symbol = n_OFDM_symbols_after_last_DRS_symbol + 1;
        end

        % construct matrix
        derotation_mat = [  0, ...                                                                                              % first symbol is not derotated
                            symbol_2_symbol_angles_normalized(1), ...                                                           % second symbol is deroteted with the angle from STF to first DRS symbols
                            repelem(symbol_2_symbol_angles_normalized(2:end), STF_and_DRS_symbol_2_symbol_spacing(2:end)), ...  % repeat the phase rotations
                            ones(1, n_OFDM_symbols_after_last_DRS_symbol)*symbol_2_symbol_angles_normalized(end)];              % symbols after the last DRS symbol

        assert(numel(derotation_mat) == N_PACKET_symb);

        % this line can be used to debug the angle
        %angles_derotation_mat = rad2deg(derotation_mat);
    
        % of course, the derotation angle from with OFDM symbols
        derotation_mat = cumsum(derotation_mat);
    
        %% create derotation matrix and derotate
        
        derotation_mat = exp(1i*(derotation_mat));
        derotation_mat = repmat(derotation_mat, N_b_DFT, 1);

        antenna_streams_mapped_rev_cfo(i) = {transmit_streams_rev_i.*derotation_mat};
    end
end
