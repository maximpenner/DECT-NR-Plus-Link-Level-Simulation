function [coarse_metric_threshold_crossing_idx] = detection(pre_fft_config, samples_antenna)
    %% precalculate parameters known at the receiver

    [n_samples_antenna, N_RX] = size(samples_antenna);

    % number of STF patterns
    L = pre_fft_config.n_pattern;

    % length of the STF with oversampling
    n_samples_STF_os = pre_fft_config.n_samples_STF_b_os;

    assert(mod(n_samples_STF_os, L) == 0);

    M = n_samples_STF_os/L;

    %% run packet detection

    % For each RX antenna, we go over all samples and calculate a metric until we cross a predefined threshold_
    % We need the sample index of this threshold crossings.
    coarse_metric_threshold_crossing_idx = zeros(N_RX, 1);

    for i=1:1:N_RX
        
        % get all samples of this particular antenna
        samples_antenna_single = samples_antenna(:,i);

        % container for detection metric
        metric = zeros(n_samples_antenna - 2*n_samples_STF_os, 1);
        
        % the step must be small enough to find every STF at every SNR
        for k = 1 : pre_fft_config.detection_threshold_step : numel(metric)

            % determine the metric at this sample index
            metric(k) = lib_rx.lib_pre_fft.coarse_metric(samples_antenna_single(k : k + n_samples_STF_os - 1), M, L, pre_fft_config.detection_minimum_power_threshold);

            % the threshold must be small enough to find every STF at every SNR
            if metric(k) >= pre_fft_config.detection_threshold_value

                % save where we have crossed the threshold, this is an absolute index
                coarse_metric_threshold_crossing_idx(i) = k;

                % abort the search for this antenna
                break;
            end
        end
    end

    % With the new cover sequence, the coarse metric is far more narrow. After detection, we jump back back a few samples.
    coarse_metric_threshold_crossing_idx = coarse_metric_threshold_crossing_idx - pre_fft_config.detection_threshold_jump_pack;

    % keep the antennas at which we found a preamble
    coarse_metric_threshold_crossing_idx = coarse_metric_threshold_crossing_idx(coarse_metric_threshold_crossing_idx > 0);

    % In a real receiver, we will start the rest of the synchronization algorithm once we hit a preamble at ANY antenna.
    % This corresponds in our case to finding the smallest index across all antennas.
    % If we didn't find any preamble at all, assume index 1, which will lead to a incorrect decoding if the STO is sufficiently large.
    if isempty(coarse_metric_threshold_crossing_idx) == true
        coarse_metric_threshold_crossing_idx = 1;
    else
        coarse_metric_threshold_crossing_idx = min(coarse_metric_threshold_crossing_idx);
    end
end
