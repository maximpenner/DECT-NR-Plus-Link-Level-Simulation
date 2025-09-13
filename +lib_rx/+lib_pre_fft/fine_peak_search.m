function [N_eff_TX_report, max_idx_fine] = fine_peak_search(verbosity, ...
                                                            pre_fft_config, ...
                                                            N_b_DFT, ...
                                                            samples_antenna, ...
                                                            samples_antenna_required, ...
                                                            STF_templates, ...
                                                            oversampling, ...
                                                            coarse_sync_idx, ...
                                                            fractional_cfo_report, ...
                                                            integer_cfo_report)
    %% precalculate parameters known at the receiver

    [~, N_RX] = size(samples_antenna);

    % number of STF patterns
    L = pre_fft_config.n_pattern;

    % length of the STF with oversampling
    n_samples_STF_os = pre_fft_config.n_samples_STF_b_os;

    assert(mod(n_samples_STF_os, L) == 0);

    %% perform crosscorrelation with STF templates for N_eff_TX and fine synchronization pointer determination

    % we have a coarse sync point, we have also determined the fractional + integer CFO

    % first create a time base which we will use to correct the fractional + integer CFO
    time_base = 0:1:(2*pre_fft_config.fine_peak_search_area + n_samples_STF_os - 1);
    time_base = time_base';

    % N_eff_TX = {1,2,4,8}, for each we will check one STF, search range can be limited by radio device class
    metric_maxima = zeros(4, N_RX);
    metric_maxima_index = zeros(4, N_RX);

    % every antenna has the same search range
    search_area_min = max(1, coarse_sync_idx - pre_fft_config.fine_peak_search_area);
    search_area_max = coarse_sync_idx + pre_fft_config.fine_peak_search_area;
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
            metric = abs(lib_util.xcorr2(samples_antenna_single_search_range, STF_template_candidate));

            % this maximum is local
            [max_val, max_idx] = max(abs(metric));

            % debugging
            if verbosity > 2
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
end
