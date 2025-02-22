function [tx, ...
          samples_antenna_tx_matlab, ...
          samples_antenna_tx_matlab_resampled, ...
          samples_antenna_tx_cpp, ...
          dect_samp_rate, ...
          hw_samp_rate] = TX_compare_numerically(ffn, json_struct)

    % these are all the information we need to generate a TX packet in Matlab and compare the IQ samples with C++ results
    mac_meta_tx.u                   = json_struct.u;
    mac_meta_tx.b                   = json_struct.b;
    mac_meta_tx.PacketLengthType    = json_struct.PacketLengthType;
    mac_meta_tx.PacketLength        = json_struct.PacketLength;
    mac_meta_tx.tm_mode_0_to_11     = json_struct.tm_mode;
    mac_meta_tx.mcs_index           = json_struct.mcs_index;
    mac_meta_tx.Z                   = json_struct.Z;
    mac_meta_tx.oversampling        = json_struct.oversampling;
    mac_meta_tx.codebook_index      = json_struct.codebook_index;
    mac_meta_tx.PLCF_type           = json_struct.PLCF_type;
    mac_meta_tx.rv                  = json_struct.rv;
    mac_meta_tx.network_id          = de2bi(json_struct.network_id ,32,'left-msb');

    iq_phase_rad                                = json_struct.tx_meta.iq_phase_rad;
    iq_phase_increment_s2s_post_resampling_rad  = json_struct.tx_meta.iq_phase_increment_s2s_post_resampling_rad;

    % read PLCF and TB bits from C++ (user bits pre channel coding)
    PLCF_bits_cpp = json_struct.data.binary.PLCF;
    TB_bits_cpp = json_struct.data.binary.TB;

    % read PCC and PDC bits from C++ (bits post channel coding)
    PCC_cpp = json_struct.data.binary.PCC;
    PDC_cpp = json_struct.data.binary.PDC;

    % matlab's function lteRateMatchTurbo() fails for BPSK modulation with large transport blocks, unknown why
    if mac_meta_tx.mcs_index == 0
        fprintf("\tMSC index is %f.\n", mac_meta_tx.mcs_index);
    end

    % create an equivalent matlab tx
    verbose = 0;
    tx = dect_tx(verbose, mac_meta_tx);

    % try to regenerate tx samples with the user bits from C++
    try
        samples_antenna_tx_matlab = tx.generate_packet(PLCF_bits_cpp, TB_bits_cpp);
    catch
        tx = [];
        samples_antenna_tx_matlab = [];
        samples_antenna_tx_matlab_resampled = [];
        samples_antenna_tx_cpp = [];
        dect_samp_rate = [];
        hw_samp_rate = [];

        fprintf("\tPacket generation failed, likely lteRateMatchTurbo() threw. Check error message. Returning.\n");
        return;
    end

    % compare PCC bit errors after coding
    err_pcc_vec = abs(PCC_cpp - double(tx.packet_data.pcc_enc_dbg.d));
    err_pcc = sum(err_pcc_vec);
    if(err_pcc > 0)
        error("File with filename %s failed at PCC comparison.", filename);
    end

    % compare PDC bit errors after coding
    err_pdc_vec = abs(PDC_cpp - double(tx.packet_data.pdc_enc_dbg.d));
    err_pdc = sum(err_pdc_vec);
    if(err_pdc > 0)
        error("File with filename %s failed at PDC comparison.", filepath_and_name);
    end

    % read tx samples from json file from C++
    samples_antenna_tx_cpp = json_struct.data.IQ.real + 1i*json_struct.data.IQ.imag;

    % json concatenates antenna samples, reshape into matrix of size N_TX x N_samples_per_antenna
    N_samples_per_antenna = numel(samples_antenna_tx_cpp)/tx.phy_4_5.tm_mode.N_TX;
    samples_antenna_tx_cpp = reshape(samples_antenna_tx_cpp, N_samples_per_antenna, tx.phy_4_5.tm_mode.N_TX);

    % resampling?
    if json_struct.resampling.L == 1 && json_struct.resampling.M == 1

        % hw sample rate and DECT sample rate are the same
        hw_samp_rate                = json_struct.resampling.samp_rate;
        dect_samp_rate              = hw_samp_rate;

        samples_antenna_tx_matlab_resampled = samples_antenna_tx_matlab;

    else

        L                           = json_struct.resampling.L;
        M                           = json_struct.resampling.M;
        hw_samp_rate                = json_struct.resampling.samp_rate;
        % note: f_pass_norm and f_stop_norm include oversampling
        f_pass_norm                 = json_struct.resampling.f_pass_norm;
        f_stop_norm                 = json_struct.resampling.f_stop_norm;
        passband_ripple_dB          = json_struct.resampling.passband_ripple_dB;
        stopband_attenuation_dB     = json_struct.resampling.stopband_attenuation_dB;
        % note: f_pass_norm and f_stop_norm include oversampling
        oversampling_minimum        = 1;

        % what is the DECT2020-NR sample rate before resampling?
        dect_samp_rate              = hw_samp_rate/L*M;

        [fir_coef_tx, ~] = lib_rx.resampling_filter_design( L, ...
                                                            M, ...
                                                            f_pass_norm, ...
                                                            f_stop_norm, ...
                                                            passband_ripple_dB, ...
                                                            stopband_attenuation_dB);

        % resample
        samples_antenna_tx_matlab_resampled = lib_rx.resampling_polyphase(samples_antenna_tx_matlab, L, M, fir_coef_tx);

        % sanity check: make sure both files are of the same length
        len_cpp = numel(samples_antenna_tx_cpp(:,1));
        len_mat = numel(samples_antenna_tx_matlab_resampled(:,1));
        if len_cpp < len_mat
            error("C++ packet shorter than matlab packet.");
        end
        len_diff = len_cpp - len_mat;
        if abs(len_cpp - len_mat) > L
            error("Packet length difference to large.");
        end

        % adjust length of matlab packet
        samples_antenna_tx_matlab_resampled = [samples_antenna_tx_matlab_resampled; zeros(len_diff, tx.phy_4_5.tm_mode.N_TX)];
    end

    % mix signal regenerated by matlab
    t_normalized = (0:1:(size(samples_antenna_tx_matlab_resampled, 1) -1))';
    mixer = exp(1i*(iq_phase_rad + iq_phase_increment_s2s_post_resampling_rad*t_normalized));
    for i=1:1:tx.phy_4_5.tm_mode.N_TX
        samples_antenna_tx_matlab_resampled(:, i) = samples_antenna_tx_matlab_resampled(:, i).*mixer;
    end

    N_samples_packet_no_GI_os_rs = json_struct.N_samples_packet_no_GI_os_rs;

    assert(0.8 <= N_samples_packet_no_GI_os_rs / numel(samples_antenna_tx_matlab_resampled(:,1)));

    % In C++, the packet is terminated at zero with the first GI samples. We only compare non-zero samples.
    A = samples_antenna_tx_matlab_resampled(1:N_samples_packet_no_GI_os_rs);
    B = samples_antenna_tx_cpp(1:N_samples_packet_no_GI_os_rs);

    % compare samples numerically
    diff = A - B;
    diff_abs = abs(diff);
    err = max(diff_abs, [], 'all');
    if(err > 0.001)

        % At this point the samples are incorrect, so there is a bug.
        % To debug, we can take the IQ samples from C++ and try to decode the frame.
        % It will likely fail, but we can still find unexpected number in the time-frequency grid.
        % Depending on the exact TX configuration, the internals of this function must be adjusted (e.g. number of RX antennas).

        figure()
        plot(abs(samples_antenna_tx_matlab_resampled))
        hold on
        plot(abs(samples_antenna_tx_cpp))
        legend('matlab', 'C++');

        error("File with filename %s failed at samples comparison. This can be due to slightly different types of mixing in C++ and Matlab.", ffn);
    end

    % power calculation
    power_per_antenna = rms(samples_antenna_tx_cpp, 1).^2;
    power_all_antennas = sum(power_per_antenna);
    fprintf("\tpower_per_antenna: ");
    for ant_idx = 1:1:numel(power_per_antenna)
        fprintf("%f ", power_per_antenna(ant_idx));
    end
    fprintf(" power_all_antennas: %f\n", power_all_antennas);

    % Define power limits for sanity check.
    % Power can be quite low when packet is short as we also include the GI.
    % Furthermore, the STF across all antennas can be low depending of the beamforming matrix.
    power_all_antennas_HIGH = 1.5;
    if numel(power_per_antenna) == 1
        power_all_antennas_LOW = 0.7;
    elseif numel(power_per_antenna) == 2
        power_all_antennas_LOW = 0.5;
    elseif numel(power_per_antenna) == 4
        power_all_antennas_LOW = 0.3;
    else
        power_all_antennas_LOW = 0.0;
        warning("Power across all antennas for 8 antennas ill-defined in ETSI standard.");
    end

    % sanity check
    if power_all_antennas < power_all_antennas_LOW || power_all_antennas > power_all_antennas_HIGH
        error("Power across all antennas out-of-bound");
    end
end