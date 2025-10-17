classdef rx_t < matlab.mixin.Copyable
    
    properties
        % copied from dectnrp_tx.tx_t()
        config;
        derived;

        % number of antennas at the receiver
        N_RX

        % processing after the FFT in frequency domain
        sto_fractional_config           % based on STF and optionally DRS
        sto_residual_config             % based on DRS
        cfo_residual_config             % based on DRS
        channel_estimation_config       % based on optionally STF and DRS
        equalization_detection_config

        packet_data;    % intermediate results during packet decoding

        HARQ_buf_40;    % PCC 40 bits
        HARQ_buf_80;    % PCC 80 bits
        HARQ_buf;       % PDC
        
        wiener;         % wiener filter for channel estimation, depends on channel environment
    end
    
    methods
        function obj = rx_t(tx)
            assert(isa(tx, "dectnrp_tx.tx_t"));
            assert(tx.config.is_valid());

            obj.config = tx.config;
            obj.derived = tx.derived;
            
            obj = obj.set_example_values();

            obj.packet_data = [];

            obj.clear_harq_buffers();
            
            obj.set_wiener(obj.channel_estimation_config.noise_estim, ...
                           obj.channel_estimation_config.f_d_hertz, ...
                           obj.channel_estimation_config.tau_rms_sec);
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
            obj.wiener = dectnrp_rx.channel_estimation_wiener_weights(obj.derived.physical_resource_mapping_DRS_cell, ...
                                                                      obj.channel_estimation_config.N_closest_DRS_pilots, ...
                                                                      obj.derived.numerology.N_b_DFT, ...
                                                                      obj.derived.N_PACKET_symb, ...
                                                                      obj.derived.numerology.N_b_CP, ...
                                                                      obj.derived.numerology.B_u_b_DFT, ...
                                                                      noise_estim, ...
                                                                      f_d_hertz, ...
                                                                      tau_rms_sec);
        end
        
        % We pass on the samples as they are received at the antennas and we try to extract the PCC and PDC bits.
        % It gets more complicated when using HARQ as calls to this function depend on each other.
        % The handle to tx is required to determine the SINR
        function [plcf_bits_recovered, tb_bits_recovered] = demod_decode_packet(obj, samples_antenna_rx)
            
            %% for the purpose of readability, extract all variables that are necessary at this stage

            verbosity           = obj.config.verbosity;

            HARQ_buf_40_        = obj.HARQ_buf_40;
            HARQ_buf_80_        = obj.HARQ_buf_80;
            HARQ_buf_           = obj.HARQ_buf;
            wiener_             = obj.wiener;
        
            mode_0_to_11        = obj.derived.tm_mode.mode_0_to_11;
            N_SS                = obj.derived.tm_mode.N_SS;
            N_eff_TX            = obj.derived.tm_mode.N_eff_TX;

            mcs                 = obj.derived.mcs;

            N_b_DFT             = obj.derived.numerology.N_b_DFT;            
            N_b_CP              = obj.derived.numerology.N_b_CP;

            N_TB_bits           = obj.derived.N_TB_bits;
            N_PACKET_symb       = obj.derived.N_PACKET_symb;
            k_b_OCC             = obj.derived.k_b_OCC;

            u                   = obj.config.u;
            b                   = obj.config.b;
            Z                   = obj.config.Z;
            network_id          = obj.config.network_id;
            PLCF_type           = obj.config.PLCF_type;
            rv                  = obj.config.rv;
            oversampling        = obj.config.oversampling;

            physical_resource_mapping_PCC_cell = obj.derived.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.derived.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.derived.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.derived.physical_resource_mapping_DRS_cell;

            assert(size(samples_antenna_rx, 1) == obj.derived.n_packet_samples*oversampling, "incorrect number of input samples");

            %% revert cover sequence by reapplying it
            samples_antenna_rx = dectnrp_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_rx, ...
                                                                                        u, ...
                                                                                        b*oversampling);

            %% OFDM demodulation a.k.a FFT
            % Switch back to frequency domain.
            % We use one version with subcarriers from oversampling removed and one with them still occupied.
            % The second version might be necessary if we have a very large, uncorrected integer CFO.
            [antenna_streams_mapped_rev, ~]= dectnrp_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion_rev(samples_antenna_rx, ...
                                                                                                                             k_b_OCC, ...
                                                                                                                             N_PACKET_symb, ...
                                                                                                                             obj.N_RX, ...
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
            if ~isempty(obj.sto_fractional_config)
                [antenna_streams_mapped_rev, sto_fractional] = dectnrp_rx.sto_fractional(antenna_streams_mapped_rev, ...
                                                                                         physical_resource_mapping_STF_cell, ...
                                                                                         obj.N_RX, ...
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
            if ~isempty(obj.cfo_residual_config)
                antenna_streams_mapped_rev = dectnrp_rx.cfo_residual(antenna_streams_mapped_rev, ...
                                                                     physical_resource_mapping_DRS_cell, ...
                                                                     physical_resource_mapping_STF_cell, ...
                                                                     obj.N_RX, ...
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
            ch_estim = dectnrp_rx.channel_estimation_wiener(antenna_streams_mapped_rev, ...
                                                            physical_resource_mapping_DRS_cell, ...
                                                            wiener_, ...
                                                            obj.N_RX, ...
                                                            N_eff_TX);
            
            %% equalization and detection
            % With the known transmission mode and the channel estimate, we can now extract the binary data from the packet.
            % Equalization here is understood as the inversion of channel effects at the subcarriers, usually paired with a symbol demapper.
            % For MIMO detection with N_SS > 1 and symbol detection, equalization is not performed explicitly but is implicit part of the underlying algorithm.
            if ~isempty(obj.equalization_detection_config)
                
                % SISO
                if N_eff_TX == 1 && obj.N_RX ==1

                    x_PCC_rev = dectnrp_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell);
                    x_PDC_rev = dectnrp_rx.equalization_SISO_zf(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell);

                % SIMO (Maximum Ratio Combining (MRC))
                elseif N_eff_TX == 1 && obj.N_RX > 1

                    x_PCC_rev = dectnrp_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PCC_cell, obj.N_RX);
                    x_PDC_rev = dectnrp_rx.equalization_SIMO_mrc(antenna_streams_mapped_rev, ch_estim, physical_resource_mapping_PDC_cell, obj.N_RX);

                % MISO (Alamouti) and MIMO (Alamouti + MRC and other modes)
                else
                    % Transmit diversity precoding
                    if ismember(mode_0_to_11, [1,5,10]) == true
                        x_PCC_rev = dectnrp_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, obj.N_RX, N_eff_TX, physical_resource_mapping_PCC_cell);
                        x_PDC_rev = dectnrp_rx.equalization_MISO_MIMO_alamouti_mrc(antenna_streams_mapped_rev, ch_estim, obj.N_RX, N_eff_TX, physical_resource_mapping_PDC_cell);
                    % MIMO modes with more than one spatial stream
                    else
                        error("MIMO modes with N_SS>1 not implemented yet.");
                    end
                end
                
            % no equalization, a great debugging tool when using AWGN channel
            else
                x_PCC_rev = dectnrp_7_transmission_encoding.subcarrier_demapping_PCC(antenna_streams_mapped_rev, physical_resource_mapping_PCC_cell);
                x_PDC_rev = dectnrp_7_transmission_encoding.subcarrier_demapping_PDC(antenna_streams_mapped_rev, physical_resource_mapping_PDC_cell);                
            end

            %% PCC and PDC decoding
            % decode PCC
            [plcf_bits_recovered, ...
             PCC_HARQ_buf_40_report, ...
             PCC_HARQ_buf_80_report, ...
             CL_report, ...
             BF_report, ...
             pcc_dec_dbg] =  dectnrp_7_transmission_encoding.PCC_decoding(x_PCC_rev, ...
                                                                          obj.config.PLCF_type, ...
                                                                          HARQ_buf_40_, ...
                                                                          HARQ_buf_80_);
            
            % decode PDC
            [tb_bits_recovered, PDC_HARQ_buf_report, pdc_dec_dbg] = dectnrp_7_transmission_encoding.PDC_decoding(x_PDC_rev, ...
                                                                                                                 N_TB_bits, ...
                                                                                                                 Z, ...
                                                                                                                 network_id, ...
                                                                                                                 PLCF_type, ...
                                                                                                                 rv, ...
                                                                                                                 mcs, ...
                                                                                                                 N_SS, ...
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
            obj.packet_data.samples_antenna_rx = samples_antenna_rx;
            obj.packet_data.ch_estim = ch_estim;
            obj.packet_data.x_PCC_rev = x_PCC_rev;
            obj.packet_data.x_PDC_rev = x_PDC_rev;
            obj.packet_data.plcf_bits_recovered = plcf_bits_recovered;
            obj.packet_data.CL_report = CL_report;
            obj.packet_data.BF_report = BF_report;
            obj.packet_data.pcc_dec_dbg = pcc_dec_dbg;
            obj.packet_data.tb_bits_recovered = tb_bits_recovered;
            obj.packet_data.pdc_dec_dbg = pdc_dec_dbg;

            if verbosity >= 1
                obj.plot_channel_estimation();
                obj.plot_scatter();
            end
        end

        % The following methods must only be called after a packet was generated by TX and processed by RX.

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
                sinr = NaN;
                return;
            end

            noise_power_measurement = obj.packet_data.x_PDC_rev - tx.packet_data.x_PDC;
            noise_power_measurement = sum(abs(noise_power_measurement).^2)/numel(noise_power_measurement);

            sinr = 10*log10(tx_power_measurement/noise_power_measurement);
        end

        function [] = plot_channel_estimation(obj)
            assert(obj.N_RX == numel(obj.packet_data.ch_estim));
        
            % extract dimensions
            [N_subc, N_symb, N_TS] =  size(obj.packet_data.ch_estim{1});
        
            figure()
            clf()

            for i=1:obj.N_RX
        
                % channel estimation of one antenna
                ch_estim_single = obj.packet_data.ch_estim{i};
        
                for j=1:N_TS
                
                    % channel estimation of one transmit stream
                    ch_estim_single_TS = ch_estim_single(:,:,j);
        
                    subplot(obj.N_RX, N_TS, (i-1)*N_TS + j);
                    hold on;
        
                    % fill subplot
                    for k=1:N_symb
        
                        % channel estimation for one OFDM symbol
                        ch_estim_single_TS_symb = ch_estim_single_TS(:, k);
        
                        % zero guards for nice plots
                        ch_estim_single_TS_symb(1:obj.derived.numerology.n_guards_bottom) = 0;
                        ch_estim_single_TS_symb(end-obj.derived.numerology.n_guards_top+1:end) = 0;
        
                        plot(-N_subc/2 : 1 : (N_subc/2-1), mag2db(abs(ch_estim_single_TS_symb)));
                    end
        
                    title_str = "RX Ch Estim: Ant " + num2str(i-1) + " TS " + num2str(j-1) + " for " + num2str(N_symb) + " Symbols";
                    title(title_str);
                    grid on
                    grid minor
                    xlim([-N_subc/2 N_subc/2])
                    ylim([-30 10])
                    xlabel("Subcarrier Index")
                    ylabel("Absolute dB")
                end
            end
        end

        function [] = plot_scatter(obj)
            h = scatterplot(obj.packet_data.x_PDC_rev);
            title("RX Scatterplot");
            set(h, 'WindowStyle', 'Docked');
        end
    end

    methods (Hidden = true)
        function obj = set_example_values(obj)
            obj.N_RX = 2;

            % can be deactivated by leaving empty
            obj.sto_fractional_config = 1;
            obj.sto_residual_config = 1;
            obj.cfo_residual_config = 1;

            obj.channel_estimation_config.N_closest_DRS_pilots = 8;
            obj.channel_estimation_config.noise_estim = 1/10^(30/10);   % 30dB SNR
            obj.channel_estimation_config.f_d_hertz = 20;               % 20Hz Doppler
            obj.channel_estimation_config.tau_rms_sec = 363e-9;         % 363ns delay spread

            % can be deactivated by leaving empty
            obj.equalization_detection_config = 1;
        end
    end
end
