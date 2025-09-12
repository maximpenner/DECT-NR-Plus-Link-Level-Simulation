function [STF_templates] = sync_STF_template(tx_config)

    %% We save STF template which are known at the receiver.

    % We save them in time domain for fine synchronization in time domaim.
    % We save them in frequency domain for correction of the integer CFO.

    % we have one stf per N_eff_TX = {1,2,4,8}
    STF_templates.time_domain = cell(4,1);
    STF_templates.freq_domain = cell(4,1);

    %% to generate the STF templates, we have to put the transmitter into a specific mode
    tx_config.PacketLengthType = 0;  % subslots
    tx_config.PacketLength = 4;      % for some tx modes (with N_eff_TX >= 4) we need at least 20 OFDM symbols

    %% one STF for each new number of effective TX antennas
    for N_eff_TX_idx=1:1:4

        % choose modes with the correct number of effective antennas
        switch N_eff_TX_idx
            case 1
                % mode with N_eff_TX = 1
                tx_config.tm_mode_0_to_11 = 0;
            case 2
                % mode with N_eff_TX = 2
                tx_config.tm_mode_0_to_11 = 2;
            case 3
                % mode with N_eff_TX = 4
                tx_config.tm_mode_0_to_11 = 6;
            case 4
                % mode with N_eff_TX = 8
                tx_config.tm_mode_0_to_11 = 11;
        end

        %% create a dummy frame

        tx_config.verbosity = 0;
        
        % create a dummy transmitter
        tx_dummy = tx_t(tx_config);
        tx_dummy.tx_config.codebook_index = 0;

        % PCC
        if tx_config.PLCF_type == 1
            PCC_user_bits = randi([0 1], 40, 1);
        elseif tx_config.PLCF_type == 2
            PCC_user_bits = randi([0 1], 80, 1);
        end

        % PDC
        N_TB_bits = tx_dummy.phy_4_5.N_TB_bits;
        PDC_user_bits = randi([0 1], N_TB_bits, 1);

        % time domain samples
        samples_for_antenna = tx_dummy.generate_packet(PCC_user_bits, PDC_user_bits);

        %% save in time domain (oversampling included)

        % how long is the dummy's preamble?
        n_STF_samples = tx_dummy.phy_4_5.n_STF_samples;

        % with oversampling STF becomes longer
        n_STF_samples = tx_config.oversampling * n_STF_samples;

        % extract STF and save in cell
        STF_templates.time_domain(N_eff_TX_idx) = {samples_for_antenna(1:n_STF_samples, 1)};

        %% save in frequency domain (oversampling excluded, is added in IFFT stage and dropped in FFT stage)

        % extract indices and values

        N_b_DFT = tx_dummy.phy_4_5.numerology.N_b_DFT;

        % make a copy or shortness
        physical_resource_mapping_STF_cell = tx_dummy.phy_4_5.physical_resource_mapping_STF_cell;

        % indices
        k_i = cell2mat(physical_resource_mapping_STF_cell(1));
        k_i_matlab = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i);

        % values
        values = cell2mat(physical_resource_mapping_STF_cell(3));

        % create empty container and fill
        STF_symbol = zeros(N_b_DFT, 1);
        STF_symbol(k_i_matlab) = values;

        STF_templates.freq_domain(N_eff_TX_idx) = {STF_symbol};
    end
end

