function [result] = TX_RX_compare_numerically(  tx, ...
                                                samples_antenna_tx_matlab_resampled, ...
                                                tx_json_struct, ...
                                                rx_json_struct)

    %% compare IQ samples

    N_RX = rx_json_struct.N_RX;

    % read tx samples from json file from C++
    samples_antenna_rx_cpp = rx_json_struct.data.IQ.real + 1i*rx_json_struct.data.IQ.imag;

    % json concatenates antenna samples, reshape into matrix of size N_TX x N_samples_per_antenna
    N_samples_per_antenna = numel(samples_antenna_rx_cpp)/N_RX;
    samples_antenna_rx_cpp = reshape(samples_antenna_rx_cpp, N_samples_per_antenna, N_RX);

    % find smaller of both sizes
    [N_samples_matlab, N_TX_matlab] = size(samples_antenna_tx_matlab_resampled);
    [N_samples_cpp, N_RX_cpp] = size(samples_antenna_rx_cpp);
    N_samples_min = min(N_samples_matlab, N_samples_cpp);

    % this comparison only works for SISO
    if N_TX_matlab ~= 1 || N_TX_matlab ~= N_RX_cpp
        error("Comparison of TX signal with RX signal only possible for SISO.");
    end

    %% determine SNR of IQ samples
    [IQ_snr_dB, IQ_err_max, IQ_err_max_idx] = get_snr_dB(samples_antenna_tx_matlab_resampled(1:N_samples_min), samples_antenna_rx_cpp(1:N_samples_min), 1.0);

    if IQ_err_max > 1000.0
        error("TX samples deviate too much from RX samples.");
    end

    fprintf("\tMaximum deviation between TX signal and RX signal is %f and idx=%d.\n", IQ_err_max, IQ_err_max_idx);
    fprintf("\tSNR between TX signal and RX signal is %f dB.\n", IQ_snr_dB);

    %% extract values right after FFT
    FFT_cpp = rx_json_struct.data.FFT.real + 1i*rx_json_struct.data.FFT.imag;

    %% extract complex derotated values of PCC
    PCC_cmplx_derotated_vec = rx_json_struct.data.PCC_cmplx_derotated_vec.real + 1i*rx_json_struct.data.PCC_cmplx_derotated_vec.imag;

    %% extract complex derotated values of PCC
    PDC_cmplx_derotated_vec = rx_json_struct.data.PDC_cmplx_derotated_vec.real + 1i*rx_json_struct.data.PDC_cmplx_derotated_vec.imag;

    %% C++ turns bit into shorts, get correct scaling factor
    PCC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QPSK;
    if tx.phy_4_5.mcs.N_bps == 1
        PDC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QPSK;
    elseif tx.phy_4_5.mcs.N_bps == 2
        PDC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QPSK;
    elseif tx.phy_4_5.mcs.N_bps == 4
        PDC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QAM16;
    elseif tx.phy_4_5.mcs.N_bps == 6
        PDC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QAM64;
    elseif tx.phy_4_5.mcs.N_bps == 8
        PDC_scaling = rx_json_struct.data.binary.SCALE_SHORT_CONV_QAM256;
    else
        error("Unknown MCS.");
    end

    %% compare PCC BER and SNR

    % first extract bits from rx
    pcc_rx_as_short = rx_json_struct.data.binary.PCC;

    % convert to +1 and -1
    pcc_rx = sign(pcc_rx_as_short);
    pcc_tx = 2*double(tx.packet_data.pcc_enc_dbg.d) - 1;

    err_bits_pcc = abs(pcc_rx - pcc_tx);
    BER_PCC = sum(err_bits_pcc) / numel(err_bits_pcc);
    [PCC_snr_dB, ~, ~] = get_snr_dB(double(pcc_rx_as_short) / PCC_scaling, pcc_tx, 1.0);
    fprintf("\tBER of PCC is %f   SNR of PCC is %f\n", BER_PCC, PCC_snr_dB);

    %% compare PDC BER and SNR
    % first extract bits from rx
    pdc_rx_as_short = rx_json_struct.data.binary.PDC;

    % convert to +1 and -1
    pdc_rx = sign(pdc_rx_as_short);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the data from c++ are actually f bits after descrambling, not f bits
    %pdc_tx = 2*double(tx.packet_data.pdc_enc_dbg.d) - 1;
    pdc_tx = 2*double(tx.packet_data.pdc_enc_dbg.f) - 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    err_bits_pdc = abs(pdc_rx - pdc_tx);
    BER_PDC = sum(err_bits_pdc) / numel(err_bits_pdc);
    [PDC_snr_dB, ~, ~] = get_snr_dB(double(pdc_rx_as_short) / PDC_scaling, pdc_tx, 1.0);
    fprintf("\tBER of PDC is %f   SNR of PDC is %f   %s\n", BER_PDC, PDC_snr_dB, tx.phy_4_5.mcs.modulation);

    %% create results

    result.samples_antenna_rx_cpp   = samples_antenna_rx_cpp;
    result.IQ_err_max_idx          = IQ_err_max_idx;

    result.IQ_snr_dB    = IQ_snr_dB;

    result.BER_PCC      = BER_PCC;
    result.PCC_snr_dB   = PCC_snr_dB;

    result.BER_PDC      = BER_PDC;
    result.PDC_snr_dB   = PDC_snr_dB;
end

function [snr_dB, err_max, err_max_idx] = get_snr_dB(sig0, sig1, ref_power)

    err = sig0 - sig1;

    [err_max, err_max_idx] = max(abs(err));

    err_power = rms(err)^2;
    snr_dB = 10*log10(ref_power/err_power);
end
