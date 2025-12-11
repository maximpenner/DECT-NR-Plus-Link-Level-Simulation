classdef sync_t < matlab.mixin.Copyable
    
    properties
        % copied from dectnrp_tx.tx_t()
        tx_config;
        derived;

        % STF templates in time and frequency domain
        n_samples_STF_b_os;
        n_samples_STF_cp_only_b_os;
        stf_templates;

        lpf_enable;

        coarse_detection_power_threshold;
        coarse_detection_step;
        coarse_detection_threshold;
        coarse_detection_jumpback;
    
        coarse_peak_search_length;
        coarse_peak_step;
        coarse_peak_threshold;
        coarse_peak_movmean;
        
        fine_peak_search_area;
    
        cfo_fractional_enable;
    
        cfo_integer_enable;
        cfo_integer_max_deviation_subcarrier_spacings;
        cfo_integer_candidate_values;

        sync_report;
    end
    
    methods
        function obj = sync_t(tx)
            assert(isa(tx, "dectnrp_tx.tx_t"));

            obj.tx_config = tx.tx_config;
            obj.derived = tx.derived;
            
            obj = obj.set_universal_values();

            obj.sync_report = [];
        end

        function [samples_antenna_rx] = synchronize(obj, samples_antenna_ch)
            %% extract size of IQ samples received
        
            [n_samples_antenna, N_RX] = size(samples_antenna_ch);

            assert(n_samples_antenna > obj.derived.n_packet_samples*obj.tx_config.oversampling, ...
                   "for synchonization more samples than the packet size must be provided");

            %% low-pass filtering 

            % If oversampling is used, we have to remove out-of-band noise, otherwise synchronization in time domain is impaired.
            % This low pass filtered signal is only used for synchronization algorithms.
            % The signal fed to the FFT still contains the out-of-band noise, since the FFT itself is an LPF.
            if obj.lpf_enable && obj.tx_config.oversampling > 1
                samples_antenna_ch_lpf = dectnrp_sync.lpf(samples_antenna_ch, obj.tx_config.oversampling);
            else
                samples_antenna_ch_lpf = samples_antenna_ch;
            end
        
            %% detection by searching for a coarse metric threshold crossing
        
            coarse_metric_threshold_crossing_idx = dectnrp_sync.coarse_detection(obj.coarse_detection_power_threshold, ...
                                                                                 obj.coarse_detection_step, ...
                                                                                 obj.coarse_detection_threshold, ...
                                                                                 obj.coarse_detection_jumpback, ...
                                                                                 obj.n_samples_STF_b_os, ...
                                                                                 obj.derived.n_STF_pattern, ...
                                                                                 samples_antenna_ch_lpf);
        
            %% after detection, extract the range of samples required for the upcoming steps
        
            A = coarse_metric_threshold_crossing_idx + 1;
            B = A + obj.coarse_peak_search_length + obj.fine_peak_search_area + obj.n_samples_STF_b_os;
        
            % if we overreach, set index to 1 which will very likely lead to a packet error if the STO is sufficiently large
            if B > n_samples_antenna
                A = 1;
                B = A + obj.coarse_peak_search_length + obj.fine_peak_search_area + obj.n_samples_STF_b_os;
            end
            
            samples_antenna_ch_lpf_required = samples_antenna_ch_lpf(A:B, :);
        
            %% search for the coarse peak starting from the coarse metric threshold crossing
        
            coarse_peak_idx = dectnrp_sync.coarse_peak_search(obj.tx_config.verbosity, ...
                                                              obj.coarse_detection_power_threshold, ...
                                                              obj.coarse_detection_threshold, ...
                                                              obj.coarse_peak_search_length, ...
                                                              obj.coarse_peak_step, ...
                                                              obj.coarse_peak_movmean, ...
                                                              obj.coarse_peak_threshold, ...
                                                              obj.n_samples_STF_b_os, ...
                                                              obj.derived.n_STF_pattern, ...
                                                              samples_antenna_ch_lpf, ...
                                                              samples_antenna_ch_lpf_required);
        
            %% at the coarse peak, determine and correct the fractional carrier frequency offset
        
            % extract samples at coarse peak
            samples_antenna_stf_at_coarse_peak = zeros(obj.n_samples_STF_b_os, N_RX);
            for i=1:1:N_RX
                samples_antenna_stf_at_coarse_peak(:,i) = samples_antenna_ch_lpf_required(coarse_peak_idx : coarse_peak_idx + obj.n_samples_STF_b_os - 1, i);
            end
        
            % use the STFs and determine a fractional CFO
            if obj.cfo_fractional_enable == true
                [samples_antenna_stf_at_coarse_peak_cfo_fractional, cfo_fractional_report] = dectnrp_sync.cfo_fractional(samples_antenna_stf_at_coarse_peak, ...
                                                                                                                         obj.n_samples_STF_b_os, ...
                                                                                                                         obj.tx_config.u);
            else
                samples_antenna_stf_at_coarse_peak_cfo_fractional = samples_antenna_stf_at_coarse_peak;
                cfo_fractional_report = 0;
            end
        
            %% at the coarse peak, determine the integer carrier frequency offset, it will be corrected later
        
            if obj.cfo_integer_enable == true
                cfo_integer_report = dectnrp_sync.cfo_integer(obj.cfo_integer_candidate_values, ...
                                                              obj.tx_config.oversampling, ...
                                                              obj.derived.numerology.N_b_DFT, ...
                                                              obj.n_samples_STF_cp_only_b_os, ...
                                                              obj.stf_templates.freq_domain, ...
                                                              samples_antenna_stf_at_coarse_peak_cfo_fractional);

            else
                cfo_integer_report = 0;
            end
        
            %% around the coarse peak, perform a cross-correlation with the STF templates
        
            [N_eff_TX_report, max_idx_fine] = dectnrp_sync.fine_peak_search(obj.tx_config.verbosity, ...
                                                                            obj.tx_config.oversampling, ...
                                                                            obj.derived.numerology.N_b_DFT, ...
                                                                            coarse_peak_idx, ...
                                                                            cfo_fractional_report, ...
                                                                            cfo_integer_report, ...
                                                                            obj.fine_peak_search_area, ...
                                                                            obj.n_samples_STF_b_os, ...
                                                                            obj.stf_templates, ...
                                                                            samples_antenna_ch_lpf, ...
                                                                            samples_antenna_ch_lpf_required);
        
            %% convert local indices to global indices
        
            max_idx_coarse = coarse_metric_threshold_crossing_idx + coarse_peak_idx;
            max_idx_fine   = coarse_metric_threshold_crossing_idx + max_idx_fine;
        
            %% extract the packet starting at the fine peak and correct the CFO for the entire packet

            n_packet_samples_os = obj.derived.n_packet_samples*obj.tx_config.oversampling;
        
            time_base = 0:1:(n_packet_samples_os - 1);
            time_base = time_base';
            
            % output
            samples_antenna_rx = zeros(n_packet_samples_os, N_RX);
        
            % final index of packet
            max_idx_fine_plus_packet = max_idx_fine + n_packet_samples_os - 1;
        
            % make sure we don't overreach
            if max_idx_fine_plus_packet > n_samples_antenna
                max_idx_fine = n_samples_antenna - n_packet_samples_os + 1;
                max_idx_fine_plus_packet = n_samples_antenna;
            end
        
            % derotate samples of each antenna
            for i=1:1:N_RX
                total_cfo = cfo_fractional_report + cfo_integer_report * 1/(obj.derived.numerology.N_b_DFT*obj.tx_config.oversampling);
                samples_antenna_rx(:,i) = samples_antenna_ch(max_idx_fine : max_idx_fine_plus_packet, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
            end

            assert(size(samples_antenna_rx, 1) == obj.derived.n_packet_samples*obj.tx_config.oversampling);
        
            %% create output report
            obj.sync_report.max_idx_coarse  = max_idx_coarse;
            obj.sync_report.fractional_cfo  = cfo_fractional_report;
            obj.sync_report.integer_cfo     = cfo_integer_report;
            obj.sync_report.N_eff_TX        = N_eff_TX_report;
            obj.sync_report.max_idx_fine    = max_idx_fine;
        end
    end

    methods (Hidden = true)
        function obj = set_universal_values(obj)
            % what is the length of the STF and how many pattern does it contain?
            switch obj.tx_config.u
                case 1
                    n_samples_STF = (64+8) * 14/9;
                    n_pattern = 7;
                case {2,4,8}
                    n_samples_STF = (64+8) * 2;
                    n_pattern = 9;
                otherwise
                    error('Unknown u.');
            end
        
            % what is the actual size of the STF with oversampling?
            obj.n_samples_STF_b_os = n_samples_STF*obj.tx_config.b*obj.tx_config.oversampling;
            obj.n_samples_STF_cp_only_b_os = obj.n_samples_STF_b_os - 64*obj.tx_config.b*obj.tx_config.oversampling;

            assert(mod(obj.n_samples_STF_b_os, obj.derived.n_STF_pattern) == 0);

            obj.stf_templates = dectnrp_rx.stf_templates(obj.tx_config);

            %% low-pass filtering to remove out-of-band noise
            
            obj.lpf_enable = true;
        
            %% STF Detection based on auto-correlation of incoming samples
            % The following values define parameters for the STF time domain synchronization, i.e. the search for where a packet starts.
            % All values are experience values, which should for work for all DECT-2020 NR packet configurations.
            % In a real receiver, some values could be optimized for the given hardware, for instance, threshold values matching ADC resolution.
            % These variables are explained and used in +dectnrp_rx.sync_STF.m.
        
            % based on the length of the STF we can estimate the length of the coarse detection metric without noise
            n_samples_coarse_metric_first_half_no_noise = obj.n_samples_STF_b_os * (n_pattern-1) / n_pattern;
        
            % power threshold: must be low enough to detect even at very low SNRs, but high enough to avoid numerical imprecision
            obj.coarse_detection_power_threshold = 0.001;
        
            % largest step is 16*b*oversampling, i.e. one STF pattern
            obj.coarse_detection_step = 8*obj.tx_config.b*obj.tx_config.oversampling;    
        
            % coarse metric is normalized between 0 and 1.0, threshold should be low enough to detect at low SNR, but not too low to avoid false alarms due to noise
            obj.coarse_detection_threshold = 0.15;
        
            % This parameter has become necessary with the newly introduced cover sequence.
            % The cover sequence has made the coarse metric very narrow, so it can happen that the metric is detected on the falling edge after the coarse peak.
            % As a countermeasure, we jump back some samples at the detection point. This way we make sure that the coarse peak search begin BEFORE the coarse peak.
            obj.coarse_detection_jumpback = obj.n_samples_STF_b_os / n_pattern * 2;
        
            %% STO Coarse Peak Search based on auto-correlation of incoming samples
        
            % search length must be long enough to definitely contain the coarse peak
            obj.coarse_peak_search_length = round(obj.coarse_detection_jumpback + 1.2*n_samples_coarse_metric_first_half_no_noise);
        
            % when using oversampling, this step can be made larger than 1
            obj.coarse_peak_step = 1;
            obj.coarse_peak_threshold = 0.15;
        
            % additional smoothing of the coarse metric
            obj.coarse_peak_movmean = [3*obj.tx_config.b*obj.tx_config.oversampling 3*obj.tx_config.b*obj.tx_config.oversampling];
        
            %% STO fine peak search around the coarse peak, based on cross-correlation with precalculated STF templates
            
            obj.fine_peak_search_area = 24*obj.tx_config.b*obj.tx_config.oversampling;
        
            %% CFO Fractional
        
            % Schmidl-Cox for an OFDM symbol with 4 repetitions instead of 2, capture range is +/- 2 subcarriers instead of +/- 1 subcarrier
        
            obj.cfo_fractional_enable = true;
        
            %% CFO Integer
        
            obj.cfo_integer_enable = true;
            
            % get the maximum physical deviation in multiples of the subcarrier spacing in use
            obj.cfo_integer_max_deviation_subcarrier_spacings = dectnrp_sync.cfo_max(obj.tx_config.u);
            
            % saw tooth of fractional CFO correction
            if obj.cfo_integer_max_deviation_subcarrier_spacings < 2
                search_lim = 0;
            elseif obj.cfo_integer_max_deviation_subcarrier_spacings < 6
                search_lim = 4;
            elseif obj.cfo_integer_max_deviation_subcarrier_spacings < 10
                search_lim = 8;
            elseif obj.cfo_integer_max_deviation_subcarrier_spacings < 14
                search_lim = 12;
            else
                error("deviation too large");
            end
            
            % we add some more possible subcarrier deviations for testing as the fractional CFO correction is not optimal
            search_lim = search_lim + 4;
            
            % search exactly within this range
            obj.cfo_integer_candidate_values = -search_lim : 4 : search_lim;
        end
    end
end
