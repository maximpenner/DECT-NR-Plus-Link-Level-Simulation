classdef dect_tx < handle
    
    properties
        verbose;        % show data during execution: 0 false, 1 only text, 2 text + plots
        mac_meta;       % data received from MAC layer
        phy_4_5;        % data from chapter 4 and 5
        
        packet_data;    % results during packet generation
    end
    
    methods
        function obj = dect_tx(verbose_arg, mac_meta_arg)
            obj.verbose = verbose_arg;
            obj.mac_meta = mac_meta_arg;
            obj.phy_4_5 = lib_util.run_chapter_4_5(verbose_arg, mac_meta_arg);

            obj.packet_data = [];
        end
        
        % We pass the bits for the PCC and PDC, and OFDM samples for each antenna are returned.
        % This function basically does all the things described in chapter 7, based on the generic procedures of chapter 6.
        function [samples_antenna_tx] = generate_packet(obj, PCC_user_bits, PDC_user_bits)
            
            %% part 3
            % for the purpose of readability, read all variables that are necessary at this stage

            verbose_        = obj.verbose;

            mode_0_to_11    = obj.phy_4_5.tm_mode.mode_0_to_11;
            N_SS            = obj.phy_4_5.tm_mode.N_SS;
            CL              = obj.phy_4_5.tm_mode.CL;
            N_TS            = obj.phy_4_5.tm_mode.N_TS;
            N_TX            = obj.phy_4_5.tm_mode.N_TX;
            N_eff_TX        = obj.phy_4_5.tm_mode.N_eff_TX;

            modulation0     = obj.phy_4_5.mcs.modulation0;

            N_b_DFT         = obj.phy_4_5.numerology.N_b_DFT;            
            N_b_CP          = obj.phy_4_5.numerology.N_b_CP;

            N_PACKET_symb   = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC         = obj.phy_4_5.k_b_OCC;
            n_STF_samples   = obj.phy_4_5.n_STF_samples;
            
            n_total_bits        = obj.phy_4_5.n_total_bits;
            n_spectrum_occupied = obj.phy_4_5.n_spectrum_occupied;

            b                   = obj.mac_meta.b;
            u                   = obj.mac_meta.u;
            Z                   = obj.mac_meta.Z;
            codebook_index      = obj.mac_meta.codebook_index;
            network_id          = obj.mac_meta.network_id;
            PLCF_type           = obj.mac_meta.PLCF_type;
            rv                  = obj.mac_meta.rv;
            oversampling        = obj.mac_meta.oversampling;

            physical_resource_mapping_PCC_cell = obj.phy_4_5.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.phy_4_5.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.phy_4_5.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.phy_4_5.physical_resource_mapping_DRS_cell;

            %% chapter 7, based on the generic procedures of chapter 6

            % The receiver needs to know if signal is beamformed or not for channel sounding purposes.
            % 7.2
            if mode_0_to_11 == 0
                precoding_identity_matrix = true;
            else
                if codebook_index == 0
                    precoding_identity_matrix = true;
                else
                    precoding_identity_matrix = false;
                end
            end
            
            % first we calculate complex samples for PCC and PDC
            % pcc_enc_dbg and pdc_enc_dbg contain debugging information            
            [x_PCC, pcc_enc_dbg] = lib_7_Transmission_modes.PCC_encoding(PCC_user_bits, CL, precoding_identity_matrix);
            [x_PDC, pdc_enc_dbg] = lib_7_Transmission_modes.PDC_encoding(PDC_user_bits,...
                                                                            n_total_bits,...
                                                                            Z,...
                                                                            network_id,...
                                                                            PLCF_type,...
                                                                            rv,...
                                                                            modulation0);
                                                                        
            % Next we map PCC and PDC to spatial streams (ss), see Table 6.3.2-1.
            % For PCC, there is only one spatial stream.
            x_PCC_ss = {x_PCC};
            if N_SS > 1
                x_PDC_ss = lib_6_generic_procedures.Spatial_Multiplexing(x_PDC, N_SS);
            else
                x_PDC_ss = {x_PDC};
            end
                                                                        
            % Transmit diversity precoding, switching to transmit streams (ts).
            % Always applied to PCC if more than one transmit stream is used.
            % Applied to PDC only if we are in mode 1, 5 or 10 according to Table 7.2-1 with N_SS = 1.
            % Otherwise, we just directly switch from spatial streams to transmit streams according to 6.3.3.1.
            if N_TS > 1
                y_PCC_ts = lib_6_generic_procedures.Transmit_diversity_precoding(x_PCC_ss, N_TS);
            else
                y_PCC_ts = x_PCC_ss;
            end
            if ismember(mode_0_to_11, [1,5,10]) == true
                y_PDC_ts = lib_6_generic_procedures.Transmit_diversity_precoding(x_PDC_ss, N_TS);
            else
                y_PDC_ts = x_PDC_ss;
            end            
           
            % We have now arrived at the transmit streams.
            % We create one matrix of size N_b_DFT x N_PACKET_symb for each transmit stream (N_TS many).
            transmit_streams = cell(N_TS,1);
            for i=1:1:N_TS
                transmit_streams(i,1) = {zeros(N_b_DFT, N_PACKET_symb)};
            end
            
            % we then map STF, DRS, PCC and PDC into those transmit stream matrices
            transmit_streams = lib_7_Transmission_modes.subcarrier_mapping_STF(transmit_streams, physical_resource_mapping_STF_cell);
            transmit_streams = lib_7_Transmission_modes.subcarrier_mapping_DRS(transmit_streams, physical_resource_mapping_DRS_cell);            
            transmit_streams = lib_7_Transmission_modes.subcarrier_mapping_PCC(transmit_streams, physical_resource_mapping_PCC_cell, y_PCC_ts);
            transmit_streams = lib_7_Transmission_modes.subcarrier_mapping_PDC(transmit_streams, physical_resource_mapping_PDC_cell, y_PDC_ts);
            
            % Beamforming (N_eff_TX many), remember N_eff_TX = N_TS
            antenna_streams_mapped = lib_6_generic_procedures.Beamforming(transmit_streams, N_TX, codebook_index);

            % switch to time domain
            samples_antenna_tx = lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion(antenna_streams_mapped,...
                                                                                                            k_b_OCC,...
                                                                                                            N_PACKET_symb,...
                                                                                                            N_TX,...
                                                                                                            N_eff_TX,...
                                                                                                            N_b_DFT,...
                                                                                                            u,...
                                                                                                            N_b_CP,...
                                                                                                            oversampling);

            % apply STF cover sequence
            samples_antenna_tx = lib_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_tx, u, b*oversampling);

            %% save packet data for debugging
            obj.packet_data.x_PCC = x_PCC;
            obj.packet_data.x_PDC = x_PDC;
            obj.packet_data.pcc_enc_dbg = pcc_enc_dbg;
            obj.packet_data.pdc_enc_dbg = pdc_enc_dbg;
            obj.packet_data.y_PCC_ts = y_PCC_ts;
            obj.packet_data.y_PDC_ts = y_PDC_ts;
            obj.packet_data.antenna_streams_mapped = antenna_streams_mapped;

            %% sanity checks
            if numel(transmit_streams) ~= N_TS
                error('We expect N_TS = %d transmit streams, but only %d streams were generated.', N_TS, numel(transmit_streams));
            end
            if numel(antenna_streams_mapped) ~= N_TX
                error('We expect N_TX = %d beamformed streams, but only %d streams were generated.', N_TX, numel(antenna_streams));
            end
            if numel(antenna_streams_mapped) ~= N_TX
                error('We expect N_TX = %d antenna streams, but only %d streams were generated.', N_TX, numel(antenna_streams_mapped));
            end
                                                                                  
            % debugging
            if verbose_ > 0
                lib_util.check_power_t_f_domain(antenna_streams_mapped, n_spectrum_occupied, samples_antenna_tx, n_STF_samples);
            end
            
            % debugging
            if verbose_ > 1
                lib_util.plot_resource_mapping_before_antenna(antenna_streams_mapped);
                lib_util.plot_STF(samples_antenna_tx, u, N_b_DFT, oversampling);
            end
        end
    end
end

