function [STF_templates] = stf_templates(config)

    assert(isa(config, "dectnrp_tx.config_t"));

    %% We save STF template which are known at the receiver.

    % We save them in time domain for fine synchronization in time domain.
    % We save them in frequency domain for correction of the integer CFO.

    % we have one STF per N_eff_TX = {1,2,4,8}
    STF_templates.time_domain = cell(4,1);
    STF_templates.freq_domain = cell(4,1);

    %% create a copy of the TX configuration and put the copy into a specific mode

    config_cpy = copy(config);

    config_cpy.PacketLengthType = 0;   % subslots
    config_cpy.PacketLength = 4;       % for some tx modes (with N_eff_TX >= 4) we need at least 20 OFDM symbols
    config_cpy.codebook_index = 0;     % no beamforming
    config_cpy.verbosity = 0;

    %% one STF for each new number of effective TX antennas
    for N_eff_TX_idx=1:1:4

        % choose modes with the correct number of effective antennas
        switch N_eff_TX_idx
            case 1
                % mode with N_eff_TX = 1
                config_cpy.tm_mode_0_to_11 = 0;
            case 2
                % mode with N_eff_TX = 2
                config_cpy.tm_mode_0_to_11 = 2;
            case 3
                % mode with N_eff_TX = 4
                config_cpy.tm_mode_0_to_11 = 6;
            case 4
                % mode with N_eff_TX = 8
                config_cpy.tm_mode_0_to_11 = 11;
        end

        % create transmitter
        tx_local = dectnrp_tx.tx_t(config_cpy);

        % create a packet
        samples_for_antenna = tx_local.generate_random_packet();

        %% save in time domain (oversampling included)

        % with oversampling STF becomes longer
        n_STF_samples_os = tx_local.derived.n_STF_samples * config_cpy.oversampling;

        % extract STF and save in cell
        STF_templates.time_domain(N_eff_TX_idx) = {samples_for_antenna(1:n_STF_samples_os, 1)};

        %% save in frequency domain (oversampling excluded, is added in IFFT stage and dropped in FFT stage)

        % readability
        N_b_DFT = tx_local.derived.numerology.N_b_DFT;
        physical_resource_mapping_STF_cell = tx_local.derived.physical_resource_mapping_STF_cell;

        % indices
        k_i = cell2mat(physical_resource_mapping_STF_cell(1));
        k_i_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, k_i);

        % values
        values = cell2mat(physical_resource_mapping_STF_cell(3));

        % create empty container and fill
        STF_symbol = zeros(N_b_DFT, 1);
        STF_symbol(k_i_matlab) = values;

        % save in cell
        STF_templates.freq_domain(N_eff_TX_idx) = {STF_symbol};
    end
end
