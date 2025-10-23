function [coarse_peak_idx] = coarse_peak_search(verbosity, ...
                                                coarse_detection_power_threshold, ...
                                                coarse_detection_threshold, ...
                                                coarse_peak_search_length, ...
                                                coarse_peak_step, ...
                                                coarse_peak_movmean, ...
                                                coarse_peak_threshold, ...
                                                n_samples_STF_b_os, ...
                                                n_STF_pattern, ...
                                                samples_antenna_ch, ...
                                                samples_antenna_ch_required)
    %% precalculate parameters known at the receiver

    [n_samples_antenna, N_RX] = size(samples_antenna_ch);

    assert(mod(n_samples_STF_b_os, n_STF_pattern) == 0);

    M = n_samples_STF_b_os/n_STF_pattern;

    %% run coarse peak search

    % We now know we most likely found an STF, now we need to find the peak of the coarse detection metric.
    % We start searching for the peak right after the index threshold_cross_idx.
    % Peaks are determined for each rx antenna individually.

    % how strong is the detected packet?
    coarse_peak_height = zeros(N_RX, 1);

    % where was the packet detected? this is now an index relative to samples_antenna_required
    coarse_peak_idx = ones(N_RX, 1);

    for i=1:1:N_RX
        
        % get samples of this particular antenna
        samples_antenna_single = samples_antenna_ch_required(:,i);

        % container for detection metric
        metric = zeros(coarse_peak_search_length, 1);

        % we continue the search with a much smaller step size
        for m = 1 : coarse_peak_step : coarse_peak_search_length

            % We found the packet too late. This can't happen in a real receiver, but here we have a limited amount of samples.
            % Nevertheless, this should barely ever happen.
            if m + n_samples_STF_b_os-1 > n_samples_antenna
                break;
            end

            % determine the metric at this sample index
            metric(m) = dectnrp_sync.coarse_metric(samples_antenna_single(m:m+n_samples_STF_b_os-1), M, n_STF_pattern, coarse_detection_power_threshold);
        end

        % warning: if sto_config.coarse_peak.step > 1, we apply the moving average across zeros!

        % apply some filtering to smooth the coarse detection
        metric_mm = movmean(metric, coarse_peak_movmean);

        % finally we search for the maximum
        [cph,cpi] = max(metric_mm);

        % save values, these values are relative to samples_antenna_required
        coarse_peak_height(i) = cph;
        coarse_peak_idx(i) = cpi;

        % debugging
        if verbosity >= 2
            figure()
            clf()
            hold on

            subplot(2,1,1)
            plot(abs(metric));
  
            yline(coarse_detection_threshold, 'r');
            text(0, coarse_detection_threshold + 0.05, 'Detection Threshold')

            xline(cpi, 'r-.');

            str = append('Coarse Synchronization Metric for antenna ', num2str(i));
            title(str);
            ylim([0 1.2]);
            grid on

            subplot(2,1,2)
            plot(abs(metric_mm));

            yline(coarse_detection_threshold, 'r');
            text(0, coarse_detection_threshold + 0.05, 'Detection Threshold')

            xline(cpi, 'r-.');
            text(cpi,cph,'Coarse Peak')

            str = append('Coarse Synchronization Metric after movmean for antenna ', num2str(i));
            title(str);
            ylim([0 1.2]);
            grid on
        end
    end

    %% coarse peak selection
    
    % we have found a peak for each rx antenna, now we have to pick the optimal peak index

    % we consider only peaks above a threshold
    idx_peak_high = coarse_peak_height >= coarse_peak_threshold;

    % remove peaks that are too low
    coarse_peak_height = coarse_peak_height(idx_peak_high);
    coarse_peak_idx    = coarse_peak_idx(idx_peak_high);

    % it can actually happen that we have no peak above the threshold, e.g. for a false alarm
    if isempty(coarse_peak_height) == false

        % sync points are peak indices weighted by the height
        normalization = coarse_peak_height/sum(coarse_peak_height);
        coarse_peak_idx = sum( coarse_peak_idx .* normalization);
        coarse_peak_idx = floor(coarse_peak_idx);

        % make sure all indexes are within range
        if coarse_peak_idx <= 0
            coarse_peak_idx = 1;
        end
        if coarse_peak_idx >= coarse_peak_search_length
            coarse_peak_idx = coarse_peak_search_length;
        end

    % if we didn't find any packets, assume the peak at 1
    else
        coarse_peak_idx = 1;
    end
end
