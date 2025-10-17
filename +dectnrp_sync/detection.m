function [coarse_metric_threshold_crossing_idx] = detection(detection_threshold_step, ...
                                                            detection_minimum_power_threshold, ...
                                                            detection_threshold_value, ...
                                                            detection_threshold_jump_pack, ...
                                                            n_samples_STF_b_os, ...
                                                            n_STF_pattern, ...
                                                            samples_antenna_ch)
    %% precalculate parameters known at the receiver

    [n_samples_antenna_ch, N_RX] = size(samples_antenna_ch);

    assert(mod(n_samples_STF_b_os, n_STF_pattern) == 0);

    M = n_samples_STF_b_os/n_STF_pattern;

    %% run packet detection

    % For each RX antenna, we go over all samples and calculate a metric until we cross a predefined threshold_
    % We need the sample index of this threshold crossings.
    coarse_metric_threshold_crossing_idx = zeros(N_RX, 1);

    for i=1:1:N_RX
        
        % get all samples of this particular antenna
        samples_antenna_single = samples_antenna_ch(:,i);

        % container for detection metric
        metric = zeros(n_samples_antenna_ch - 2*n_samples_STF_b_os, 1);
        
        % the step must be small enough to find every STF at every SNR
        for k = 1 : detection_threshold_step : numel(metric)

            % determine the metric at this sample index
            metric(k) = dectnrp_sync.coarse_metric(samples_antenna_single(k : k + n_samples_STF_b_os - 1), M, n_STF_pattern, detection_minimum_power_threshold);

            % the threshold must be small enough to find every STF at every SNR
            if metric(k) >= detection_threshold_value

                % save where we have crossed the threshold, this is an absolute index
                coarse_metric_threshold_crossing_idx(i) = k;

                % abort the search for this antenna
                break;
            end
        end
    end

    % With the new cover sequence, the coarse metric is far more narrow. After detection, we jump back back a few samples.
    coarse_metric_threshold_crossing_idx = coarse_metric_threshold_crossing_idx - detection_threshold_jump_pack;

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
