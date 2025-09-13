function [samples_antenna_sync, sync_report] = sync(verbosity, ...
                                                    pre_fft_config, ...
                                                    u, ...
                                                    N_b_DFT, ...
                                                    samples_antenna, ...
                                                    STF_templates, ...
                                                    n_packet_samples, ...
                                                    oversampling)
    %% precalculate parameters known at the receiver

    [n_samples_antenna, N_RX] = size(samples_antenna);

    % number of STF patterns
    L = pre_fft_config.n_pattern;

    % length of the STF with oversampling
    n_samples_STF_os = pre_fft_config.n_samples_STF_b_os;

    assert(mod(n_samples_STF_os, L) == 0);

    %% packet detection by searching for a coarse metric threshold crossing

    coarse_metric_threshold_crossing_idx = lib_rx.lib_pre_fft.detection(pre_fft_config, samples_antenna);

    %% after detection, extract the range of samples required for the upcoming steps

    A = coarse_metric_threshold_crossing_idx + 1;
    B = A + pre_fft_config.coarse_peak_search_length + pre_fft_config.fine_peak_search_area + n_samples_STF_os;

    % if we overreach, set index to 1 which will very likely lead to a packet error if the STO is sufficiently large
    if B > n_samples_antenna
        A = 1;
        B = A + pre_fft_config.coarse_peak_search_length + pre_fft_config.fine_peak_search_area + n_samples_STF_os;
    end
    
    samples_antenna_required = samples_antenna(A:B, :);

    %% starting from the coarse metric threshold crossing, search for the coarse peak

    coarse_peak_idx = lib_rx.lib_pre_fft.coarse_peak_search(verbosity, ...
                                                            pre_fft_config, ...
                                                            samples_antenna, ...
                                                            samples_antenna_required);

    %% at the coarse peak, determine and correct the fractional carrier frequency offset

    % extract samples at coarse peak
    samples_antenna_stf_at_coarse_peak = zeros(n_samples_STF_os, N_RX);
    for i=1:1:N_RX
        samples_antenna_stf_at_coarse_peak(:,i) = samples_antenna_required(coarse_peak_idx : coarse_peak_idx + n_samples_STF_os - 1, i);
    end

    % use the STFs and determine a fractional CFO
    if pre_fft_config.active_fractional == true
        [samples_antenna_stf_at_coarse_peak_cfo_fractional, cfo_fractional_report] = lib_rx.lib_pre_fft.cfo_fractional(samples_antenna_stf_at_coarse_peak, ...
                                                                                                                       n_samples_STF_os, ...
                                                                                                                       u);
    else
        samples_antenna_stf_at_coarse_peak_cfo_fractional = samples_antenna_stf_at_coarse_peak;
        cfo_fractional_report = 0;
    end

    %% at the coarse peak, determine the integer carrier frequency offset, it will be corrected later

    if pre_fft_config.active_integer == true
        integer_cfo_report = lib_rx.lib_pre_fft.cfo_integer(pre_fft_config, ...
                                                            samples_antenna_stf_at_coarse_peak_cfo_fractional, ...
                                                            STF_templates.freq_domain, ...
                                                            N_b_DFT, ...
                                                            oversampling);
    else
        integer_cfo_report = 0;
    end

    %% around the coarse peak, perform a cross-correlation with the STF templates

    [N_eff_TX_report, max_idx_fine] = lib_rx.lib_pre_fft.fine_peak_search(verbosity, ...
                                                                          pre_fft_config, ...
                                                                          N_b_DFT, ...
                                                                          samples_antenna, ...
                                                                          samples_antenna_required, ...
                                                                          STF_templates, ...
                                                                          oversampling, ...
                                                                          coarse_peak_idx, ...
                                                                          cfo_fractional_report, ...
                                                                          integer_cfo_report);

    %% convert local indices to global indices

    max_idx_coarse = coarse_metric_threshold_crossing_idx + coarse_peak_idx;
    max_idx_fine   = coarse_metric_threshold_crossing_idx + max_idx_fine;

    %% extract the packet starting at the fine peak and correct the CFO for the entire packet
    
    % packet is longer when oversampled
    n_packet_samples_os = n_packet_samples*oversampling;

    % we need a new time base with the length of one frame
    time_base = 0:1:(n_packet_samples_os - 1);
    time_base = time_base';
    
    % output container
    samples_antenna_sync = zeros(n_packet_samples_os, N_RX);

    % final index of packet
    max_idx_fine_plus_packet = max_idx_fine + n_packet_samples_os - 1;

    % make sure we don't overreach
    if max_idx_fine_plus_packet > n_samples_antenna
        max_idx_fine = n_samples_antenna - n_packet_samples_os + 1;
        max_idx_fine_plus_packet = n_samples_antenna;
    end

    % derotate samples of each antenna
    for i=1:1:N_RX
        total_cfo = cfo_fractional_report + integer_cfo_report * 1/(N_b_DFT*oversampling);
        samples_antenna_sync(:,i) = samples_antenna(max_idx_fine : max_idx_fine_plus_packet, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
    end

    %% create output report
    sync_report.max_idx_coarse   = max_idx_coarse;
    sync_report.fractional_cfo   = cfo_fractional_report;
    sync_report.integer_cfo      = integer_cfo_report;
    sync_report.N_eff_TX         = N_eff_TX_report;
    sync_report.max_idx_fine     = max_idx_fine;
end
