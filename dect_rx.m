classdef dect_rx < handle
    
    properties
        verbose;        % show data during execution: 0 false, 1 only text, 2 text + plots
        mac_meta;       % data received from MAC layer
        phy_4_5;        % data from chapter 4 and 5

        packet_data;    % intermediate results during packet decoding
        
        tx_handle;      % handle to tx for debugging
        ch_handle;      % handle to rf channel for debugging

        STF_templates;  % for synchronization based on STF, both time and frequency domain
        harq_buf_40;    % these are the cbsbuffers used by matlab for harq operation, PCC 40 bits
        harq_buf_80;    % these are the cbsbuffers used by matlab for harq operation, PCC 80 bits
        harq_buf;       % these are the cbsbuffers used by matlab for harq operation, PDC
        
        wiener;         % structure with data about wiener filter for channel estimation, depends on channel environment
    end
    
    methods
        function obj = dect_rx(verbose_arg, mac_meta_arg)
            % if channel estimation type is not specifically defined, use Wiener as a default
            if ~isfield(mac_meta_arg, 'active_ch_estim_type')
                mac_meta_arg.active_ch_estim_type = 'wiener';
            end

            % if equalization/detection is not specifically turned off, activate if by default
            if ~isfield(mac_meta_arg, 'active_equalization_detection')
                mac_meta_arg.active_equalization_detection = true;
            end

            obj.verbose = verbose_arg;
            obj.mac_meta = mac_meta_arg;
            obj.phy_4_5 = lib_util.run_chapter_4_5(verbose_arg, mac_meta_arg);

            obj.packet_data = [];

            obj.tx_handle = [];
            obj.ch_handle = [];

            % we only load STF templates if they are needed
            if mac_meta_arg.synchronization.pre_FFT.active == true
                obj.STF_templates = lib_rx.sync_STF_template(mac_meta_arg);
            end

            obj.harq_buf_40 = [];
            obj.harq_buf_80 = [];
            obj.harq_buf = [];
            
            % init wiener coefficients with reasonable default values
            obj.overwrite_wiener(1/10^(30/10), ...  % 30dB SNR
                                 20, ...            % 20Hz Doppler
                                 363e-9);           % 363ns delay spread
        end

        % Wiener filter interpolates between DRS pilots. Optimal interpolation depends on channel properties.
        % We use a static filter which is used regardless of the instantaneous channel.
        % For noise, the best case value is assumed, for delay and Doppler spread the worst case value.
        % To improve performance, different sets should be precalculated for different SNRs.
        function [] = overwrite_wiener(obj, noise_estim, f_d_hertz, tau_rms_sec)
            % We estimate the channel for each subcarrier based on the closest DRS pilots in time and frequency.
            % This value defines how many of the closest DRS pilots we consider. Ideally, this value should depend on the SNR. 
            N_closest_DRS_pilots = 8;

            obj.wiener = lib_rx.channel_estimation_wiener_weights(  obj.phy_4_5.physical_resource_mapping_DRS_cell,...
                                                                    N_closest_DRS_pilots,...
                                                                    obj.phy_4_5.numerology.N_b_DFT,...
                                                                    obj.phy_4_5.N_PACKET_symb,...
                                                                    obj.phy_4_5.numerology.N_b_CP, ...
                                                                    obj.phy_4_5.numerology.B_u_b_DFT,...
                                                                    noise_estim, ...
                                                                    f_d_hertz, ...
                                                                    tau_rms_sec);
        end
        
        % We pass on the samples as they are received at the antennas and we try to extract the PCC and PDC bits.
        % It gets more complicated when using HARQ as calls to this function depend on each other.
        function [PCC_user_bits, PDC_user_bits] = demod_decode_packet(obj, samples_antenna_rx)

            % for the purpose of readability, extract all variables that are necessary at this stage

            verbose_            = obj.verbose;
            tx_handle_          = obj.tx_handle;
            STF_templates_      = obj.STF_templates;
            harq_buf_40_        = obj.harq_buf_40;
            harq_buf_80_        = obj.harq_buf_80;
            harq_buf_           = obj.harq_buf;
            wiener_             = obj.wiener;
        
            mode_0_to_11        = obj.phy_4_5.tm_mode.mode_0_to_11;
            N_eff_TX            = obj.phy_4_5.tm_mode.N_eff_TX;

            modulation0         = obj.phy_4_5.mcs.modulation0;

            N_b_DFT             = obj.phy_4_5.numerology.N_b_DFT;            
            N_b_CP              = obj.phy_4_5.numerology.N_b_CP;

            N_TB_bits           = obj.phy_4_5.N_TB_bits;
            N_PACKET_symb       = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC             = obj.phy_4_5.k_b_OCC;

            n_packet_samples    = obj.phy_4_5.n_packet_samples;

            b                   = obj.mac_meta.b;
            u                   = obj.mac_meta.u;
            Z                   = obj.mac_meta.Z;
            network_id          = obj.mac_meta.network_id;
            PLCF_type           = obj.mac_meta.PLCF_type;
            rv                  = obj.mac_meta.rv;
            N_RX                = obj.mac_meta.N_RX;
            oversampling        = obj.mac_meta.oversampling;

            synchronization     = obj.mac_meta.synchronization;

            active_ch_estim_type = obj.mac_meta.active_ch_estim_type;

            active_equalization_detection = obj.mac_meta.active_equalization_detection;

            physical_resource_mapping_PCC_cell = obj.phy_4_5.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.phy_4_5.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.phy_4_5.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.phy_4_5.physical_resource_mapping_DRS_cell;

            %% sync of integer STO and fractional + integer CFO in time domain before FFT, and extraction of N_eff_TX into STO_CFO_report
            if synchronization.pre_FFT.active == true

                % The property snr_db of the RF channel refers to the in-band noise.
                % If oversampling is used, we have to remove out-of-band noise, otherwise synchronization in time domain is impaired.
                if oversampling > 1
                    % Kaiser LPF
                    lowpassfilter = designfilt( 'lowpassfir', ...
                                                'PassbandFrequency', 0.6 / oversampling, ...
                                                'StopbandFrequency',0.8 / oversampling, ...
                                                'PassbandRipple', 10, ...
                                                'StopbandAttenuation', 30, ...
                                                'SampleRate', 1, ...
                                                'DesignMethod', 'kaiserwin', ...
                                                'MinOrder', 'even');
                    
                    assert(mod(numel(lowpassfilter.Coefficients), 2) == 1, 'fractional LPF delay');
                    
                    samples_antenna_rx = filter(lowpassfilter, samples_antenna_rx);
                    
                    % compensate deterministic filter delay prior to synchronization
                    lowpassfilter_delay = (numel(lowpassfilter.Coefficients)-1)/2;
                    samples_antenna_rx(1:end-lowpassfilter_delay) = samples_antenna_rx(lowpassfilter_delay+1:end);
                end

                % The number of samples received is larger than the number of samples in a packet.
                % In this function, we try to synchronize the packet and extract the exact number of samples in the packet, which we assume to be known.
                % In a real receiver, the number of samples is unknown. We would first have to decode the PCC which lies in the first few OFDM symbols and extract that information.
                [samples_antenna_rx_sto_cfo, STO_CFO_report] = lib_rx.sync_STF( verbose_,...
                                                                                u,...
                                                                                N_b_DFT,...
                                                                                samples_antenna_rx,...
                                                                                STF_templates_,...
                                                                                n_packet_samples,...
                                                                                oversampling,...
                                                                                synchronization.pre_FFT.sto_config,...
                                                                                synchronization.pre_FFT.cfo_config);

                obj.packet_data.STO_CFO_report = STO_CFO_report;
            else
                % assume the input samples are synchronized and have the correct length
                samples_antenna_rx_sto_cfo = samples_antenna_rx;

                % add empty report
                obj.packet_data.STO_CFO_report = [];
            end

            %% revert cover sequence by reapplying it
            samples_antenna_rx_sto_cfo = lib_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_rx_sto_cfo, u, b*oversampling);

            %% OFDM demodulation a.k.a FFT
            % Switch back to frequency domain.
            % We use one version with subcarriers from oversampling removed and one with them still occupied.
            % The second version might be necessary if we have a very large, uncorrected integer CFO.
            [antenna_streams_mapped_rev, ~]= lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion_rev(   samples_antenna_rx_sto_cfo,...
                                                                                                                            k_b_OCC,...
                                                                                                                            N_PACKET_symb,...
                                                                                                                            N_RX,...
                                                                                                                            N_eff_TX,...
                                                                                                                            N_b_DFT,...
                                                                                                                            u,...
                                                                                                                            N_b_CP,...
                                                                                                                            oversampling);

            %% remove fractional + residual STO based on STF and DRS
            % Idea: Delay in time domain leads to increasing phase rotation within an OFDM symbol, but steady across packet.
            % This is known as the phase error gradient (PEG).
            % Residual STO may also be due to a Symbol Clock Offset (SCO).
            % ToDo: add DRS to fractional STO estimation (optional, increases latency in a real receiver)
            % ToDo: residual STO based on STF and DRS
            if synchronization.post_FFT.sto_fractional == true
                [antenna_streams_mapped_rev, sto_fractional] = lib_rx.sync_STO_fractional(antenna_streams_mapped_rev, physical_resource_mapping_STF_cell, N_RX, oversampling);

                % add to report
                obj.packet_data.STO_CFO_report.sto_fractional = sto_fractional;
            else
                % add empty
                obj.packet_data.STO_CFO_report.sto_fractional = 0;
            end

            %% remove residual CFO correction based on STF and DRS
            % Idea: CFO leads to steady phase rotation within an OFDM symbol, but increasing phase rotation across packet.
            % This is known as the common phase error (CPE).
            % Is a real receiver, a CPU can also be caused by phase noise.
            % ToDo: add STF to residual CFO estimation
            if synchronization.post_FFT.cfo_residual == true
                antenna_streams_mapped_rev = lib_rx.sync_CFO_residual(  antenna_streams_mapped_rev,...
                                                                        physical_resource_mapping_DRS_cell,...
                                                                        physical_resource_mapping_STF_cell,...
                                                                        N_RX,...
                                                                        N_eff_TX);
            end

            %% PCC decoding and packet data extraction
            % Decode PCC and extract all relevant data about the packet.
            % We know N_eff_TX from STO sync, so we also know the DRS pattern.
            % With STF and DRS known, we can now estimate the channel for PCC, which lies at the beginning of each packet.
            % PCC is always transmitted with transmit diversity coding.
            % For now it is done down below together with the PDC.
            % ToDo
                                                                                                                    
            %% channel estimation
            % Now that we know the length of the packet from the PCC, we can determine a channel estimate.
            % Output is a cell(N_RX,1), each cell with a matrix of size N_b_DFT x N_PACKET_symb x N_TX.
            % For subcarriers which are unused the channel estimate can be NaN or +/- infinity.
            if strcmp(active_ch_estim_type,'wiener') == true

                % real world channel estimation based on precalculated wiener filter coefficients assuming worst-case channel conditions
                ch_estim = lib_rx.channel_estimation_wiener(antenna_streams_mapped_rev, physical_resource_mapping_DRS_cell, wiener_, N_RX, N_eff_TX);

            elseif strcmp(active_ch_estim_type,'least squares') == true

                % Zero forcing at the pilots and linear interpolation for subcarriers inbetween.
                % It works without channel knowledge, but has a fairly large error floor.
                % The error floor is comes from only considering the closest neighbours.
                ch_estim = lib_rx.channel_estimation_ls(antenna_streams_mapped_rev, physical_resource_mapping_STF_cell, physical_resource_mapping_DRS_cell, N_RX, N_eff_TX);


            else
                error('Unknown channel estimation type %s.', active_ch_estim_type);
            end
            
            %% equalization and detection
            % With the known transmission mode and the channel estimate, we can now extract the binary data from the packet.
            % Equalization here is understood as the inversion of channel effects at the subcarriers, usually paired with a symbol demapper.
            % For MIMO detection with N_SS > 1 and symbol detection, equalization is not performed explicitly but is implicit part of the underlying algorithm.
            if active_equalization_detection == true
                
                % SISO
                if N_eff_TX == 1 && N_RX ==1

                    x_PCC_rev = lib_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell);
                    x_PDC_rev = lib_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell);

                % SIMO (Maximum Ratio Combining (MRC))
                elseif N_eff_TX == 1 && N_RX > 1

                    x_PCC_rev = lib_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell, N_RX);
                    x_PDC_rev = lib_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell, N_RX);

                % MISO (Alamouti) and MIMO (Alamouti + MRC and other modes)
                else
                    % Transmit diversity precoding
                    if ismember(mode_0_to_11, [1,5,10]) == true
                        x_PCC_rev = lib_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, N_RX, N_eff_TX, physical_resource_mapping_PCC_cell);
                        x_PDC_rev = lib_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, N_RX, N_eff_TX, physical_resource_mapping_PDC_cell);
                    % MIMO modes with more than one spatial stream
                    else
                        error("MIMO modes with N_SS>1 not implemented yet.");
                    end
                end
                
            % no equalization, a great debugging tool when using awgn channel
            else
                x_PCC_rev = lib_7_Transmission_modes.subcarrier_unmapping_PCC(antenna_streams_mapped_rev, physical_resource_mapping_PCC_cell);
                x_PDC_rev = lib_7_Transmission_modes.subcarrier_unmapping_PDC(antenna_streams_mapped_rev, physical_resource_mapping_PDC_cell);                
            end

            %% PCC and PDC decoding
            % decode PCC
            [PCC_user_bits, PCC_harq_buf_40_report,...
                            PCC_harq_buf_80_report,...
                            CL_report,...
                            BF_report,...
                            PCC_report,...
                            pcc_dec_dbg] =  lib_7_Transmission_modes.PCC_decoding(  x_PCC_rev, ...
                                                                                    harq_buf_40_,...
                                                                                    harq_buf_80_);
            
            % decode PDC
            [PDC_user_bits, PDC_harq_buf_report, pdc_dec_dbg] = lib_7_Transmission_modes.PDC_decoding(  x_PDC_rev,...
                                                                                                        N_TB_bits,...
                                                                                                        Z,...
                                                                                                        network_id,...
                                                                                                        PLCF_type,...
                                                                                                        rv,...
                                                                                                        modulation0,...
                                                                                                        harq_buf_);

            %% HARQ buffers
            % save harq buffer for next call
            if numel(PCC_harq_buf_40_report) > 0
                obj.harq_buf_40 = PCC_harq_buf_40_report;
            end
            if numel(PCC_harq_buf_80_report) > 0
                obj.harq_buf_80 = PCC_harq_buf_80_report;
            end            
            if numel(PDC_harq_buf_report) > 0
                obj.harq_buf = PDC_harq_buf_report;
            end

            %% save packet data for debugging
            obj.packet_data.CL_report   = CL_report;
            obj.packet_data.BF_report   = BF_report;
            obj.packet_data.PCC_report  = PCC_report;
            obj.packet_data.pcc_dec_dbg = pcc_dec_dbg;
            obj.packet_data.pdc_dec_dbg = pdc_dec_dbg;

            % debugging
            if verbose_ > 0
                disp('##### RX debugging ######');
                
                % compare complex symbols of PCC and PDF
                fprintf('Max error distance PCC: %.20f\n', max(abs(x_PCC_rev - tx_handle_.packet_data.x_PCC)));
                fprintf('Max error distance PDC: %.20f\n', max(abs(x_PDC_rev - tx_handle_.packet_data.x_PDC)));

                % measure the sinr in dB by comparing the complex PDC symbols
                tx_power_measurement = tx_handle_.packet_data.x_PDC;
                tx_power_measurement = sum(abs(tx_power_measurement).^2)/numel(tx_power_measurement);
                noise_power_measurement = x_PDC_rev - tx_handle_.packet_data.x_PDC;
                noise_power_measurement = sum(abs(noise_power_measurement).^2)/numel(noise_power_measurement);
                sinr = 10*log10(tx_power_measurement/noise_power_measurement);
                fprintf('Measured f domain SINR: %f dB\n', sinr);
                fprintf('Measured sto_fractional: %f samples\n', obj.packet_data.STO_CFO_report.sto_fractional);

                if verbose_ > 1
                    lib_dbg.plot_chestim(ch_estim, obj.phy_4_5.numerology.n_guards_bottom, obj.phy_4_5.numerology.n_guards_top);
                    scatterplot(x_PDC_rev);                    
                end
            end
        end
    end
end

