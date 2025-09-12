classdef rx_t < matlab.mixin.Copyable
    
    properties
        tx_config;      % data received from MAC layer
        phy_4_5;        % data from chapter 4 and 5
        packet_data;    % intermediate results during packet decoding

        rx_config;      % additional configuration of the receiver

        HARQ_buf_40;    % PCC 40 bits
        HARQ_buf_80;    % PCC 80 bits
        HARQ_buf;       % PDC
        
        wiener;         % wiener filter for channel estimation, depends on channel environment
    end
    
    methods
        function obj = rx_t(tx, rx_config)
            assert(isa(tx, "lib_types.tx_t"));
            assert(isa(rx_config, "lib_types.rx_config_t"));

            % create copies
            obj.tx_config = tx.tx_config;
            obj.phy_4_5 = tx.phy_4_5;
            obj.packet_data = [];

            obj.rx_config = rx_config;

            obj.clear_harq_buffers();
            
            obj.set_wiener(obj.rx_config.channel_estimation_config.noise_estim, ...
                           obj.rx_config.channel_estimation_config.f_d_hertz, ...
                           obj.rx_config.channel_estimation_config.tau_rms_sec);
        end

        function [] = clear_harq_buffers(obj)
            obj.HARQ_buf_40 = [];
            obj.HARQ_buf_80 = [];
            obj.HARQ_buf = [];
        end

        % Wiener filter interpolates between DRS pilots. Optimal interpolation depends on channel properties.
        % We use a static filter which is used regardless of the instantaneous channel.
        % For noise, the best case value is assumed, for delay and Doppler spread the worst case value.
        % To improve performance, different sets should be precalculated for different SNRs.
        function [] = set_wiener(obj, noise_estim, f_d_hertz, tau_rms_sec)
            obj.wiener = lib_rx.channel_estimation_wiener_weights(obj.phy_4_5.physical_resource_mapping_DRS_cell, ...
                                                                  obj.rx_config.channel_estimation_config.N_closest_DRS_pilots, ...
                                                                  obj.phy_4_5.numerology.N_b_DFT, ...
                                                                  obj.phy_4_5.N_PACKET_symb, ...
                                                                  obj.phy_4_5.numerology.N_b_CP, ...
                                                                  obj.phy_4_5.numerology.B_u_b_DFT, ...
                                                                  noise_estim, ...
                                                                  f_d_hertz, ...
                                                                  tau_rms_sec);
        end
        
        % We pass on the samples as they are received at the antennas and we try to extract the PCC and PDC bits.
        % It gets more complicated when using HARQ as calls to this function depend on each other.
        % The handle to tx is required to determine the SINR
        function [plcf_bits_recovered, tb_bits_recovered] = demod_decode_packet(obj, samples_antenna_rx)
            %% for the purpose of readability, extract all variables that are necessary at this stage

            verbosity           = obj.tx_config.verbosity;

            HARQ_buf_40_        = obj.HARQ_buf_40;
            HARQ_buf_80_        = obj.HARQ_buf_80;
            HARQ_buf_           = obj.HARQ_buf;
            wiener_             = obj.wiener;

            stf_templates       = obj.rx_config.stf_templates;
        
            mode_0_to_11        = obj.phy_4_5.tm_mode.mode_0_to_11;
            N_eff_TX            = obj.phy_4_5.tm_mode.N_eff_TX;

            modulation0         = obj.phy_4_5.mcs.modulation0;

            N_b_DFT             = obj.phy_4_5.numerology.N_b_DFT;            
            N_b_CP              = obj.phy_4_5.numerology.N_b_CP;

            N_TB_bits           = obj.phy_4_5.N_TB_bits;
            N_PACKET_symb       = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC             = obj.phy_4_5.k_b_OCC;

            n_packet_samples    = obj.phy_4_5.n_packet_samples;

            u                   = obj.tx_config.u;
            b                   = obj.tx_config.b;
            Z                   = obj.tx_config.Z;
            network_id          = obj.tx_config.network_id;
            PLCF_type           = obj.tx_config.PLCF_type;
            rv                  = obj.tx_config.rv;
            N_RX                = obj.rx_config.N_RX;
            oversampling        = obj.tx_config.oversampling;

            pre_fft_config      = obj.rx_config.pre_fft_config;

            equalization_detection_config = obj.rx_config.equalization_detection_config;

            physical_resource_mapping_PCC_cell = obj.phy_4_5.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.phy_4_5.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.phy_4_5.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.phy_4_5.physical_resource_mapping_DRS_cell;

            %% synchronization before the FFT based on the STF
            if ~isempty(pre_fft_config)

                assert(size(samples_antenna_rx, 1) > obj.phy_4_5.n_packet_samples*oversampling, ...
                       "for synchonization more samples than the packet size must be provided");

                % If oversampling is used, we have to remove out-of-band noise.
                % Otherwise synchronization in time domain is impaired.
                if oversampling > 1
                    samples_antenna_rx = lib_rx.lib_pre_fft.lpf(samples_antenna_rx, oversampling);
                end

                % The number of samples received is larger than the number of samples in a packet.
                % We try to synchronize the packet and extract the exact number of samples in the packet, which we assume to be known.
                % In a real receiver, the number of samples is unknown.
                % We would first have to decode the PCC which lies in the first few OFDM symbols and extract that information.
                [samples_antenna_rx_sto_cfo, sync_report] = lib_rx.lib_pre_fft.sync(verbosity, ...
                                                                                    pre_fft_config, ...
                                                                                    u, ...
                                                                                    N_b_DFT, ...
                                                                                    samples_antenna_rx, ...
                                                                                    stf_templates, ...
                                                                                    n_packet_samples, ...
                                                                                    oversampling);

                obj.packet_data.sync_report = sync_report;
            else
                assert(size(samples_antenna_rx, 1) == obj.phy_4_5.n_packet_samples*oversampling, ...
                       "without synchonization the same number of samples as the packet contains must be provided");

                % assume the input samples are synchronized and have the correct length
                samples_antenna_rx_sto_cfo = samples_antenna_rx;

                obj.packet_data.sync_report = [];
            end

            %% revert cover sequence by reapplying it
            samples_antenna_rx_sto_cfo = lib_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_rx_sto_cfo, ...
                                                                                            u, ...
                                                                                            b*oversampling);

            %% OFDM demodulation a.k.a FFT
            % Switch back to frequency domain.
            % We use one version with subcarriers from oversampling removed and one with them still occupied.
            % The second version might be necessary if we have a very large, uncorrected integer CFO.
            [antenna_streams_mapped_rev, ~]= lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion_rev(samples_antenna_rx_sto_cfo, ...
                                                                                                                         k_b_OCC, ...
                                                                                                                         N_PACKET_symb, ...
                                                                                                                         N_RX, ...
                                                                                                                         N_eff_TX, ...
                                                                                                                         N_b_DFT, ...
                                                                                                                         u, ...
                                                                                                                         N_b_CP, ...
                                                                                                                         oversampling);

            %% remove fractional + residual STO based on STF and DRS
            % Idea: A delay in time domain leads to increasing phase rotation within an OFDM symbol, but steady across the packet.
            % This is known as the phase error gradient (PEG).
            % A residual STO may also be due to a Symbol Clock Offset (SCO).
            % ToDo: residual STO based on STF plus DRS
            if ~isempty(obj.rx_config.post_fft_sto_fractional_config)
                [antenna_streams_mapped_rev, sto_fractional] = lib_rx.sto_fractional(antenna_streams_mapped_rev, ...
                                                                                     physical_resource_mapping_STF_cell, ...
                                                                                     N_RX, ...
                                                                                     oversampling);

                % add to report
                obj.packet_data.residual_report.sto_fractional = sto_fractional;
            else
                % add empty
                obj.packet_data.residual_report.sto_fractional = 0;
            end

            %% remove residual CFO correction based on STF and DRS
            % Idea: A CFO leads to steady phase rotation within an OFDM symbol, but increasing phase rotation across the packet.
            % This is known as the common phase error (CPE).
            % Is a real receiver, a CFO can also be caused by phase noise.
            % ToDo: residual CFO based on STF plus DRS
            if ~isempty(obj.rx_config.post_FFT_cfo_residual_config)
                antenna_streams_mapped_rev = lib_rx.cfo_residual(antenna_streams_mapped_rev, ...
                                                                 physical_resource_mapping_DRS_cell, ...
                                                                 physical_resource_mapping_STF_cell, ...
                                                                 N_RX, ...
                                                                 N_eff_TX);
            end

            %% PCC decoding and packet data extraction
            % Decode PCC and extract all relevant data about the packet.
            % We know N_eff_TX from STO sync, so we also know the DRS pattern.
            % With STF and DRS known, we can now estimate the channel for PCC, which lies at the beginning of each packet.
            % PCC is always transmitted with transmit diversity coding.
            % For now it is done down below together with the PDC.
                                                                                                                    
            %% channel estimation
            % Now that we know the length of the packet from the PCC, we can determine a channel estimate.
            % Output is a cell(N_RX,1), each cell with a matrix of size N_b_DFT x N_PACKET_symb x N_TX.
            % For subcarriers which are unused the channel estimate can be NaN or +/- infinity.
            ch_estim = lib_rx.channel_estimation_wiener(antenna_streams_mapped_rev, ...
                                                        physical_resource_mapping_DRS_cell, ...
                                                        wiener_, ...
                                                        N_RX, ...
                                                        N_eff_TX);
            
            %% equalization and detection
            % With the known transmission mode and the channel estimate, we can now extract the binary data from the packet.
            % Equalization here is understood as the inversion of channel effects at the subcarriers, usually paired with a symbol demapper.
            % For MIMO detection with N_SS > 1 and symbol detection, equalization is not performed explicitly but is implicit part of the underlying algorithm.
            if ~isempty(equalization_detection_config)
                
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
                x_PCC_rev = lib_7_transmission_encoding.subcarrier_demapping_PCC(antenna_streams_mapped_rev, physical_resource_mapping_PCC_cell);
                x_PDC_rev = lib_7_transmission_encoding.subcarrier_demapping_PDC(antenna_streams_mapped_rev, physical_resource_mapping_PDC_cell);                
            end

            %% PCC and PDC decoding
            % decode PCC
            [plcf_bits_recovered, ...
             PCC_HARQ_buf_40_report, ...
             PCC_HARQ_buf_80_report, ...
             CL_report, ...
             BF_report, ...
             pcc_dec_dbg] =  lib_7_transmission_encoding.PCC_decoding(x_PCC_rev, ...
                                                                      obj.tx_config.PLCF_type, ...
                                                                      HARQ_buf_40_, ...
                                                                      HARQ_buf_80_);
            
            % decode PDC
            [tb_bits_recovered, PDC_HARQ_buf_report, pdc_dec_dbg] = lib_7_transmission_encoding.PDC_decoding(x_PDC_rev, ...
                                                                                                             N_TB_bits, ...
                                                                                                             Z, ...
                                                                                                             network_id, ...
                                                                                                             PLCF_type, ...
                                                                                                             rv, ...
                                                                                                             modulation0, ...
                                                                                                             HARQ_buf_);

            %% save HARQ buffers for next call
            if numel(PCC_HARQ_buf_40_report) > 0
                obj.HARQ_buf_40 = PCC_HARQ_buf_40_report;
            end
            if numel(PCC_HARQ_buf_80_report) > 0
                obj.HARQ_buf_80 = PCC_HARQ_buf_80_report;
            end            
            if numel(PDC_HARQ_buf_report) > 0
                obj.HARQ_buf = PDC_HARQ_buf_report;
            end

            %% save packet data
            obj.packet_data.x_PCC_rev = x_PCC_rev;
            obj.packet_data.x_PDC_rev = x_PDC_rev;
            obj.packet_data.CL_report = CL_report;
            obj.packet_data.BF_report = BF_report;
            obj.packet_data.plcf_bits_recovered = plcf_bits_recovered;
            obj.packet_data.tb_bits_recovered = tb_bits_recovered;
            obj.packet_data.pcc_dec_dbg = pcc_dec_dbg;
            obj.packet_data.pdc_dec_dbg = pdc_dec_dbg;

            if verbosity > 0
                disp('##### RX Packet ######');
                fprintf('Measured sto_fractional: %f samples\n\n', obj.packet_data.residual_report.sto_fractional);
            end

            if verbosity > 1
                lib_util.lib_plot_print.plot_channel_estimation(ch_estim, ...
                                                                obj.phy_4_5.numerology.n_guards_bottom, ...
                                                                obj.phy_4_5.numerology.n_guards_top);
                scatterplot(x_PDC_rev);                    
            end
        end

        function n = get_pcc_bit_errors_uncoded(obj, tx)
            n = sum(abs(double(tx.packet_data.pcc_enc_dbg.d) - double(obj.packet_data.pcc_dec_dbg.d_hard)));
        end

        function n = get_pdc_bit_errors_uncoded(obj, tx)
            n = sum(abs(double(tx.packet_data.pdc_enc_dbg.d) - double(obj.packet_data.pdc_dec_dbg.d_hard)));
        end

        function ret = are_plcf_bits_equal(obj, tx)
            % in case decoding failed, i.e. incorrect CRC
            if numel(obj.packet_data.plcf_bits_recovered) == 0
                ret = false;
                return;
            end

            % decoding was successful, however, the CRC can still be incorrect
            if sum(tx.packet_data.plcf_bits - double(obj.packet_data.plcf_bits_recovered)) == 0
                ret = true;
            else
                ret = false;
            end
        end

        function ret = are_tb_bits_equal(obj, tx)
            % in case decoding failed, i.e. incorrect CRC
            if numel(obj.packet_data.tb_bits_recovered) == 0
                ret = false;
                return;
            end

            % decoding was successful, however, the CRC can still be incorrect
            if sum(tx.packet_data.tb_bits - double(obj.packet_data.tb_bits_recovered)) == 0
                ret = true;
            else
                ret = false;
            end
        end

        function sinr = get_sinr(obj, tx)
            tx_power_measurement = tx.packet_data.x_PDC;
            tx_power_measurement = sum(abs(tx_power_measurement).^2)/numel(tx_power_measurement);

            if numel(obj.packet_data.x_PDC_rev) == 0
                sinr = NaN
                return;
            end

            noise_power_measurement = obj.packet_data.x_PDC_rev - tx.packet_data.x_PDC;
            noise_power_measurement = sum(abs(noise_power_measurement).^2)/numel(noise_power_measurement);

            sinr = 10*log10(tx_power_measurement/noise_power_measurement);
        end
    end
end
