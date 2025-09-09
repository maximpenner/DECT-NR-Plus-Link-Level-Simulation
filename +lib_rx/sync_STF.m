function [samples_antenna_sto_cfo, STO_CFO_report] = sync_STF(  verbose,...
                                                                u,...
                                                                N_b_DFT,...
                                                                samples_antenna,...
                                                                STF_templates,...
                                                                n_packet_samples,...
                                                                oversampling,...
                                                                sto_config,...
                                                                cfo_config)

    %% precalculate parameters known at the receiver a-priori

    % packet size
    [n_samples_ch, N_RX] = size(samples_antenna);

    % STF pattern
    L = sto_config.n_pattern;

    % what is the length of the STF with oversampling included?
    n_samples_STF_os = sto_config.n_samples_STF_b_os;

    % sanity check
    if mod(n_samples_STF_os, L) ~= 0
        error("STF length not a multiple of L.");
    end

    % what is the length of a single pattern?
    M = n_samples_STF_os/L;

    %% detection by coarse metric threshold crossing

    % For each RX antenna, we go over all samples and calculate a metric until we cross a predefined threshold.
    % We need the sample index of this threshold crossings.
    coarse_threshold_crossing_idx = zeros(N_RX, 1);

    for i=1:1:N_RX
        
        % get all samples of this particular antenna
        samples_antenna_single = samples_antenna(:,i);

        % container for detection metric
        metric = zeros(n_samples_ch - 2*n_samples_STF_os, 1);
        
        % the step must be small enough to find every STF at every SNR
        for k = 1 : sto_config.threshold.step : numel(metric)

            % determine the metric at this sample index
            metric(k) = lib_rx.sync_STO_coarse_metric(samples_antenna_single(k : k + n_samples_STF_os - 1), M, L, sto_config.minimum_power_threshold);

            % the threshold must be small enough to find every STF at every SNR
            if metric(k) >= sto_config.threshold.value

                % save where we have crossed the threshold, this is an absolute index
                coarse_threshold_crossing_idx(i) = k;

                % abort the search for this antenna
                break;
            end
        end
    end

    % With the new cover sequence, the coarse metric is far more narrow. After detection, we jump back back a few samples to make sure we find
    % the coarse peak.
    coarse_threshold_crossing_idx = coarse_threshold_crossing_idx - sto_config.threshold.jump_pack;

    % keep the antennas at which we found a preamble
    coarse_threshold_crossing_idx = coarse_threshold_crossing_idx(coarse_threshold_crossing_idx > 0);

    % In a real receiver, we will start the rest of the synchronization algorithm once we hit a preamble at ANY antenna.
    % This corresponds in our case to finding the smallest index across all antennas.
    % If we didn't find any preamble at all, assume index 1, which will lead to a incorrect decoding if the STO is sufficiently large.
    if isempty(coarse_threshold_crossing_idx) == true
        coarse_threshold_crossing_idx = 1;
    else
        coarse_threshold_crossing_idx = min(coarse_threshold_crossing_idx);
    end

    %% after detection, extract the range of samples what will be required for the upcoming synchronization steps

    req_start = coarse_threshold_crossing_idx + 1;
    req_end = req_start + sto_config.coarse_peak.search_length + sto_config.fine.search_area + n_samples_STF_os;

    % if we overreach, set index to 1 which will very likely to an erroneous packet (if the STO is sufficiently large)
    if req_end > n_samples_ch
        req_start = 1;
        req_end = req_start + sto_config.coarse_peak.search_length + sto_config.fine.search_area + n_samples_STF_os;
    end
    
    samples_antenna_required = samples_antenna(req_start:req_end, :);

    %% coarse peak search

    % We now know we most likely found an STF, now we need to find the peak of the coarse detection metric.
    % We start searching for the peak right after the index threshold_cross_idx.
    % Peaks are determined for each rx antenna individually.

    % how strong is the detected packet?
    coarse_peak_height = zeros(N_RX, 1);

    % where was the packet detected? this is now an index relative to samples_antenna_required
    coarse_peak_idx = ones(N_RX, 1);

    for i=1:1:N_RX
        
        % get samples of this particular antenna
        samples_antenna_single = samples_antenna_required(:,i);

        % container for detection metric
        metric = zeros(sto_config.coarse_peak.search_length, 1);

        % we continue the search with a much smaller step size
        for m = 1 : sto_config.coarse_peak.step : sto_config.coarse_peak.search_length

            % We found the packet too late. This can't happen in a real receiver, but here we have a limited amount of samples.
            % Nevertheless, this should barely ever happen.
            if m + n_samples_STF_os-1 > n_samples_ch
                break;
            end

            % determine the metric at this sample index
            metric(m) = lib_rx.sync_STO_coarse_metric(samples_antenna_single(m:m+n_samples_STF_os-1), M, L, sto_config.minimum_power_threshold);
        end

        % warning: if sto_config.coarse_peak.step > 1, we apply the moving average across zeros!

        % apply some filtering to smooth the coarse detection
        metric_mm = movmean(metric, sto_config.coarse_peak.movmean);

        % finally we search for the maximum
        [cph,cpi] = max(metric_mm);

        % save values, these values are relative to samples_antenna_required
        coarse_peak_height(i) = cph;
        coarse_peak_idx(i) = cpi;

        % debugging
        if verbose > 1
            figure()
            clf()
            hold on

            subplot(2,1,1)
            plot(abs(metric));
  
            yline(sto_config.threshold.value, 'r');
            text(0,sto_config.threshold.value + 0.05,'Detection Threshold')

            xline(cpi, 'r-.');

            str = append('Coarse Synchronization Metric for antenna ', num2str(i));
            title(str);
            ylim([0 1.2]);
            grid on

            subplot(2,1,2)
            plot(abs(metric_mm));

            yline(sto_config.threshold.value, 'r');
            text(0,sto_config.threshold.value + 0.05,'Detection Threshold')

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
    idx_peak_high = coarse_peak_height >= sto_config.coarse_peak.threshold;

    % remove peaks that are too low
    coarse_peak_height = coarse_peak_height(idx_peak_high);
    coarse_peak_idx    = coarse_peak_idx(idx_peak_high);

    % it can actually happen that we have no peak above the threshold, e.g. for a false alarm
    if isempty(coarse_peak_height) == false

        % sync points are peak indices weighted by the height
        normalization = coarse_peak_height/sum(coarse_peak_height);
        coarse_sync_idx = sum( coarse_peak_idx .* normalization);
        coarse_sync_idx = floor(coarse_sync_idx);

        % make sure all indexes are within range
        if coarse_sync_idx <= 0
            coarse_sync_idx = 1;
        end
        if coarse_sync_idx >= sto_config.coarse_peak.search_length
            coarse_sync_idx = sto_config.coarse_peak.search_length;
        end

    % if we didn't find any packets, assume synchronization at 1
    else
        coarse_sync_idx = 1;
    end

    %% fractional cfo calculation

    % create a container with the coarsly synchronized STFs only
    samples_antenna_stf_at_coarse_peak = zeros(n_samples_STF_os, N_RX);
    for i=1:1:N_RX
        samples_antenna_stf_at_coarse_peak(:,i) = samples_antenna_required(coarse_sync_idx : coarse_sync_idx + n_samples_STF_os - 1, i);
    end

    if cfo_config.active_fractional == true
        % use the STFs and determine a fractional CFO
        [samples_antenna_stf_at_coarse_peak_post_fractional,...
            fractional_cfo_report] = lib_rx.sync_CFO_fractional(verbose, samples_antenna_stf_at_coarse_peak, n_samples_STF_os, u, cfo_config);
    else
        samples_antenna_stf_at_coarse_peak_post_fractional = samples_antenna_stf_at_coarse_peak;
        fractional_cfo_report = 0;
    end

    %% integer cfo calculation

    if cfo_config.active_integer == true
        integer_cfo_report = lib_rx.sync_CFO_integer(   samples_antenna_stf_at_coarse_peak_post_fractional,...
                                                        STF_templates.freq_domain,...
                                                        N_b_DFT,...
                                                        oversampling,...
                                                        sto_config,...
                                                        cfo_config);
    else
        integer_cfo_report = 0;
    end

    %% perform crosscorrelation with STF templates for N_eff_TX and fine synchronization pointer determination

    % we have a coarse sync point, we have also determined the fractional + integer CFO

    % first create a time base which we will use to correct the fractional + integer CFO
    time_base = 0:1:(2*sto_config.fine.search_area + n_samples_STF_os - 1);
    time_base = time_base';

    % N_eff_TX = {1,2,4,8}, for each we will check one STF, search range can be limited by radio device class
    metric_maxima = zeros(4, N_RX);
    metric_maxima_index = zeros(4, N_RX);

    % every antenna has the same search range
    search_area_min = max(1, coarse_sync_idx - sto_config.fine.search_area);
    search_area_max = coarse_sync_idx + sto_config.fine.search_area;
    n_search_area = search_area_max - search_area_min + n_samples_STF_os;

    for i=1:1:N_RX
        
        % get the samples of this particular antenna
        samples_antenna_single = samples_antenna_required(:,i);
        
        % try every possible stf type
        for j=1:1:numel(STF_templates.time_domain)

            % stf templates in time domain are already oversampled
            STF_template_candidate = cell2mat(STF_templates.time_domain(j));

            % normalize power of stf
            STF_template_candidate = STF_template_candidate/rms(STF_template_candidate);

            % extract the search range in which we will be looking for the fine sync peaks
            samples_antenna_single_search_range = samples_antenna_single(search_area_min : search_area_max + n_samples_STF_os - 1, :);

            % derotate samples in search range
            total_cfo = fractional_cfo_report + integer_cfo_report * 1/(N_b_DFT*oversampling);
            samples_antenna_single_search_range = samples_antenna_single_search_range.*exp(1i*2*pi*(-total_cfo)*time_base(1:n_search_area));

            % perform a cross correlation between the samples and the stf template
            metric = abs(lib_rx.sync_xcorr(samples_antenna_single_search_range, STF_template_candidate));

            % this maximum is local
            [max_val, max_idx] = max(abs(metric));

            % debugging
            if verbose > 2
                figure()

                plot(abs(metric));

                xline(max_idx, 'r-.');
                text(max_idx,max_val,'Fine Peak')

                str = append('Fine Synchronization Metric for antenna ', num2str(i), ' and STF ', num2str(j));
                title(str);
                grid on
            end

            % save the maximum
            metric_maxima(j,i) = abs(metric(max_idx));

            % this index is converted to global
            max_idx = max_idx + search_area_min - 1;

            % save global maximum
            metric_maxima_index(j,i) = max_idx;
        end
    end
    
    %% find best fitting N_eff_TX

    % average peaks across antennas
    metric_maxima_sum = sum(metric_maxima, 2);
    
    % find the best fitting stf
    [~, j_best] = max(metric_maxima_sum);
    
    % how many effective antennas do we have?
    switch j_best
        case 1
            N_eff_TX_report = 1;
        case 2
            N_eff_TX_report = 2;
        case 3
            N_eff_TX_report = 4;
        case 4
            N_eff_TX_report = 8;
    end

    %% fine best fitting fine synchronization point

    metric_maxima_fine = metric_maxima(j_best, :);

    max_idx_fine = metric_maxima_index(j_best, :);

    % weight a detection point by the height of the coarse metric of the respective antenna
    if sum(metric_maxima_fine) > 0
        normalization = metric_maxima_fine/sum(metric_maxima_fine);
    else
        normalization = 1/numel(metric_maxima_fine);
    end
    max_idx_fine = sum(max_idx_fine .* normalization);
    max_idx_fine = floor(max_idx_fine);

    % minimum is 1
    if max_idx_fine <= 0
        max_idx_fine = 1;
    end

    %% convert local indices to global indices

    max_idx_coarse = coarse_threshold_crossing_idx + coarse_sync_idx;
    max_idx_fine   = coarse_threshold_crossing_idx + max_idx_fine;

    %% extract packet starting at the fine synchronization point and correct the CFO for the entire packet
    
    % packet is longer when oversampled
    n_packet_samples = n_packet_samples*oversampling;

    % we need a new time base with the length of one frame
    time_base = 0:1:(n_packet_samples - 1);
    time_base = time_base';
    
    % output container
    samples_antenna_sto_cfo = zeros(n_packet_samples, N_RX);

    % final index of packet
    max_idx_fine_plus_packet = max_idx_fine + n_packet_samples - 1;

    % make sure we don't overreach
    if max_idx_fine_plus_packet > n_samples_ch
        max_idx_fine = n_samples_ch - n_packet_samples + 1;
        max_idx_fine_plus_packet = n_samples_ch;
    end

    for i=1:1:N_RX
        % derotate samples
        total_cfo = fractional_cfo_report + integer_cfo_report * 1/(N_b_DFT*oversampling);
        samples_antenna_sto_cfo(:,i) = samples_antenna(max_idx_fine : max_idx_fine_plus_packet, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
    end

    %% create output report
    STO_CFO_report.max_idx_coarse   = max_idx_coarse;
    STO_CFO_report.fractional_cfo   = fractional_cfo_report;
    STO_CFO_report.integer_cfo      = integer_cfo_report;
    STO_CFO_report.N_eff_TX         = N_eff_TX_report;
    STO_CFO_report.max_idx_fine     = max_idx_fine;
end

