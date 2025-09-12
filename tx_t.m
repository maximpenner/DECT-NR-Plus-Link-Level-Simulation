classdef tx_t < matlab.mixin.Copyable
    
    properties
        tx_config;      % data received from MAC layer
        phy_4_5;        % data from chapter 4 and 5
        packet_data;    % intermediate results during packet decoding
    end
    
    methods
        function obj = tx_t(tx_config)
            assert(isa(tx_config, "tx_config_t"));
            assert(tx_config.is_valid());

            obj.tx_config = tx_config;
            obj.phy_4_5 = lib_util.run_chapter_4_5(tx_config);
            obj.packet_data = [];

            if obj.tx_config.verbosity > 1
                lib_util.lib_plot_print.plot_resource_mapping(obj.phy_4_5.numerology, ...
                                                              obj.phy_4_5.N_PACKET_symb, ...
                                                              obj.phy_4_5.tm_mode.N_TS, ...
                                                              obj.phy_4_5.tm_mode.N_SS, ...
                                                              obj.phy_4_5.physical_resource_mapping_STF_cell, ...
                                                              obj.phy_4_5.physical_resource_mapping_DRS_cell, ...
                                                              obj.phy_4_5.physical_resource_mapping_PCC_cell, ...
                                                              obj.phy_4_5.physical_resource_mapping_PDC_cell);
            end
        end
        
        function [samples_antenna_tx] = generate_packet(obj, plcf_bits, tb_bits)
            
            %% for the purpose of readability, read all variables that are necessary at this stage

            verbosity_          = obj.tx_config.verbosity;

            mode_0_to_11        = obj.phy_4_5.tm_mode.mode_0_to_11;
            N_SS                = obj.phy_4_5.tm_mode.N_SS;
            CL                  = obj.phy_4_5.tm_mode.CL;
            N_TS                = obj.phy_4_5.tm_mode.N_TS;
            N_TX                = obj.phy_4_5.tm_mode.N_TX;
            N_eff_TX            = obj.phy_4_5.tm_mode.N_eff_TX;

            modulation0         = obj.phy_4_5.mcs.modulation0;

            N_b_DFT             = obj.phy_4_5.numerology.N_b_DFT;
            N_b_CP              = obj.phy_4_5.numerology.N_b_CP;

            N_PACKET_symb       = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC             = obj.phy_4_5.k_b_OCC;
            n_STF_samples       = obj.phy_4_5.n_STF_samples;
            
            n_total_bits        = obj.phy_4_5.n_total_bits;
            n_spectrum_occupied = obj.phy_4_5.n_spectrum_occupied;

            u                   = obj.tx_config.u;
            b                   = obj.tx_config.b;
            Z                   = obj.tx_config.Z;
            codebook_index      = obj.tx_config.codebook_index;
            network_id          = obj.tx_config.network_id;
            PLCF_type           = obj.tx_config.PLCF_type;
            rv                  = obj.tx_config.rv;
            oversampling        = obj.tx_config.oversampling;

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
            
            % First, we calculate complex samples for PCC and PDC.
            % pcc_enc_dbg and pdc_enc_dbg contain debugging information.
            [x_PCC, pcc_enc_dbg] = lib_7_transmission_encoding.PCC_encoding(plcf_bits, ...
                                                                            CL, ...
                                                                            precoding_identity_matrix);
            [x_PDC, pdc_enc_dbg] = lib_7_transmission_encoding.PDC_encoding(tb_bits, ...
                                                                            n_total_bits, ...
                                                                            Z, ...
                                                                            network_id, ...
                                                                            PLCF_type, ...
                                                                            rv, ...
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
            transmit_streams = lib_7_transmission_encoding.subcarrier_mapping_STF(transmit_streams, physical_resource_mapping_STF_cell);
            transmit_streams = lib_7_transmission_encoding.subcarrier_mapping_DRS(transmit_streams, physical_resource_mapping_DRS_cell);            
            transmit_streams = lib_7_transmission_encoding.subcarrier_mapping_PCC(transmit_streams, physical_resource_mapping_PCC_cell, y_PCC_ts);
            transmit_streams = lib_7_transmission_encoding.subcarrier_mapping_PDC(transmit_streams, physical_resource_mapping_PDC_cell, y_PDC_ts);
            
            % Beamforming (N_eff_TX many), remember N_eff_TX = N_TS
            antenna_streams_mapped = lib_6_generic_procedures.Beamforming(transmit_streams, N_TX, codebook_index);

            % switch to time domain
            samples_antenna_tx = lib_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion(antenna_streams_mapped, ...
                                                                                                         k_b_OCC, ...
                                                                                                         N_PACKET_symb, ...
                                                                                                         N_TX, ...
                                                                                                         N_eff_TX, ...
                                                                                                         N_b_DFT, ...
                                                                                                         u, ...
                                                                                                         N_b_CP, ...
                                                                                                         oversampling);

            % apply STF cover sequence
            samples_antenna_tx = lib_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_tx, u, b*oversampling);

            assert(numel(transmit_streams) == N_TS);
            assert(numel(antenna_streams_mapped) == N_TX);
            assert(numel(antenna_streams_mapped) == N_TX);       

            %% save packet data
            obj.packet_data.plcf_bits = plcf_bits;
            obj.packet_data.tb_bits = tb_bits;
            obj.packet_data.x_PCC = x_PCC;
            obj.packet_data.x_PDC = x_PDC;
            obj.packet_data.pcc_enc_dbg = pcc_enc_dbg;
            obj.packet_data.pdc_enc_dbg = pdc_enc_dbg;
            obj.packet_data.y_PCC_ts = y_PCC_ts;
            obj.packet_data.y_PDC_ts = y_PDC_ts;
            obj.packet_data.antenna_streams_mapped = antenna_streams_mapped;

            if verbosity_ > 0
                lib_util.lib_plot_print.print_power(antenna_streams_mapped, n_spectrum_occupied, samples_antenna_tx, n_STF_samples);
            end   
        end

        function [samples_antenna_tx] = generate_random_packet(obj)

            % PLCF
            if obj.tx_config.PLCF_type == 1
                plcf_bits = randi([0 1], 40, 1);
            elseif obj.tx_config.PLCF_type == 2
                plcf_bits = randi([0 1], 80, 1);
            else
                error("undefined PLCF type");
            end
            
            % transport block
            tb_bits = randi([0 1], obj.phy_4_5.N_TB_bits, 1);

            samples_antenna_tx = obj.generate_packet(plcf_bits, tb_bits);
        end
    end
end
