function [pre_fft_config] = config(u, b, oversampling)

    % what is the length of the STF in samples for u=1 and b=1 and how many pattern does it contain?
    switch u
        case 1
            n_samples_STF = (64+8) * 14/9;
            n_pattern = 7;
        case {2,4,8}
            n_samples_STF = (64+8) * 2;
            n_pattern = 9;
        otherwise
            error('Unknown u.');
    end

    % save
    pre_fft_config.n_samples_STF = n_samples_STF;
    pre_fft_config.n_pattern = n_pattern;

    % what is the actual size of the STF with oversampling?
    pre_fft_config.n_samples_STF_b_os = n_samples_STF*b*oversampling;
    pre_fft_config.n_samples_STF_cp_only_b_os = pre_fft_config.n_samples_STF_b_os - 64*b*oversampling;

    %% STF Detection based on auto-correlation of incoming samples

    % based on the length of the STF we can estimate the length of the coarse detection metric without noise
    n_samples_coarse_metric_first_half_no_noise = n_samples_STF * b * oversampling * (n_pattern-1) / n_pattern;

    % The following values define parameters for the STF time domain synchronization, i.e. the search for where a packet starts.
    % All values are experience values, which should for work for all DECT-2020 NR packet configurations.
    % In a real receiver, some values could be optimized for the given hardware, for instance, threshold values matching ADC resolution.
    % These variables are explained and used in +lib_rx.sync_STF.m.

    % power threshold: must be low enough to detect even at very low SNRs, but high enough to avoid numerical imprecision
    pre_fft_config.detection_minimum_power_threshold = 0.001;

    % largest step is 16*b*oversampling, i.e. one STF pattern
    pre_fft_config.detection_threshold_step = 8*b*oversampling;    

    % coarse metric is normalized between 0 and 1.0, threshold should be low enough to detect at low SNR, but not too low to avoid false alarms due to noise
    pre_fft_config.detection_threshold_value = 0.15;

    % This parameter has become necessary with the newly introduced cover sequence.
    % The cover sequence has made the coarse metric very narrow, so it can happen that the metric is detected on the falling edge after the coarse peak.
    % As a countermeasure, we jump back some samples at the detection point. This way we make sure that the coarse peak search begin BEFORE the coarse peak.
    pre_fft_config.detection_threshold_jump_pack = pre_fft_config.n_samples_STF_b_os / n_pattern * 2;

    %% STO Coarse Peak Search based on auto-correlation of incoming samples

    % starting from the detection point, the search length must be long enough to definitely contain the coarse peak
    pre_fft_config.coarse_peak_search_length = round(pre_fft_config.detection_threshold_jump_pack + 1.2*n_samples_coarse_metric_first_half_no_noise);

    % when using oversampling, this step can be made larger than 1
    pre_fft_config.coarse_peak_step = 1;
    pre_fft_config.coarse_peak_threshold = 0.15;

    % additional smoothing of the coarse metric
    pre_fft_config.coarse_peak_movmean = [3*b*oversampling 3*b*oversampling];

    %% STO fine peak search around the coarse peak, based on cross-correlation with precalculated STF templates
    
    pre_fft_config.fine_peak_search_area = 24*b*oversampling;

    %% CFO Fractional

    % Schmidl-Cox for an OFDM symbol with 4 repetitions instead of 2, capture range is +/- 2 subcarriers instead of +/- 1 subcarrier

    pre_fft_config.cfo_fractional_enable = true;

    %% CFO Integer
    %
    % How large is the search space for the integer CFO?
    %
    %   According to the DECT-2020 NR standard, battery powered devices are allowed to have up to 30ppm.
    %   At the receiver, we see a maximum of 2*30ppm = 60ppm deviation.
    %   At 6GHz, 60 ppm corresponds to 360kHz.
    %   The minimum subcarrier spacing is 27kHz, so the maximum CFO is 360kHz/27kHz = 13.33333 subcarrier spacing.
    %
    %   With larger u, the subcarrier spacing is larger as well and the actual possible deviation becomes smaller.

    pre_fft_config.cfo_integer_enable = true;
    
    % get the maximum physical deviation in multiples of the subcarrier spacing in use
    pre_fft_config.cfo_integer_max_deviation_subcarrier_spacings = sync_CFO_max_estimation(u);
    
    % saw tooth of fractional CFO correction
    if pre_fft_config.cfo_integer_max_deviation_subcarrier_spacings < 2
        search_lim = 0;
    elseif pre_fft_config.cfo_integer_max_deviation_subcarrier_spacings < 6
        search_lim = 4;
    elseif pre_fft_config.cfo_integer_max_deviation_subcarrier_spacings < 10
        search_lim = 8;
    elseif pre_fft_config.cfo_integer_max_deviation_subcarrier_spacings < 14
        search_lim = 12;
    end
    
    % we add some more possible subcarrier deviations for testing as the fractional CFO correction is not optimal
    search_lim = search_lim + 4;
    
    % search exactly within this range
    pre_fft_config.cfo_integer_candidate_values = -search_lim : 4 : search_lim;
end

function [cfo_integer_max_deviation_subcarrier_spacings] = sync_CFO_max_estimation(u)

    assert(ismember(u, [1 2 4 8]));

    % subcarrier spacing
    subc_spacing = 27000 * u;
    
    % maximum operational frequency for DECT-2020 NR
    f_max = 6e9;

    % according to table 5.5.2-1 in part 2, battery powered devices can have up to 30 ppm in extreme conditions
    ppm_max = 30;
    
    % number of subcarriers that our signal can deviate (2 as we can have 30 ppm at the transmitter and 30 ppm at the receiver)
    cfo_integer_max_deviation_subcarrier_spacings = f_max * 2 * ppm_max/1e6 / subc_spacing;
end
