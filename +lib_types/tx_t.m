classdef tx_t < matlab.mixin.Copyable
    
    properties
        tx_config;      % data received from MAC layer
        phy_4_5;        % data from clause 4 and 5
        packet_data;    % intermediate results during packet decoding
    end
    
    methods
        function obj = tx_t(tx_config)
            assert(isa(tx_config, "lib_types.tx_config_t"));
            assert(tx_config.is_valid());

            obj.tx_config = tx_config;
            obj.run_clause_4_5();
            obj.packet_data = [];
        end

        % This method is called in the constructor. It basically does all the calculations of clause 4 and 5.
        function [] = run_clause_4_5(obj)
            %% for the purpose of readability
            u                   = obj.tx_config.u;
            b                   = obj.tx_config.b;
            PacketLengthType    = obj.tx_config.PacketLengthType;
            PacketLength        = obj.tx_config.PacketLength;
            tm_mode_0_to_11     = obj.tx_config.tm_mode_0_to_11;
            mcs_index           = obj.tx_config.mcs_index;
            Z                   = obj.tx_config.Z;
        
            %% start generating the packet structure with the functions provided in the technical specification (TS)
        
            % 7.2
            tm_mode = lib_7_transmission_encoding.transmission_modes(tm_mode_0_to_11);
        
            % Annex A
            mcs = lib_Annex_A.modulation_and_coding_scheme(mcs_index);
        
            % 4.3
            numerology = lib_4_physical_layer_principles.numerologies(u,b);
        
            % 4.4
            [T_frame, N_FRAME_slot, T_slot] = lib_4_physical_layer_principles.frame_structure();
        
            % 4.5
            k_b_OCC = lib_4_physical_layer_principles.physical_resources(numerology.N_b_OCC);
        
            % 5.1
            N_PACKET_symb = lib_5_physical_layer_transmissions.Transmission_packet_structure(numerology, ...
                                                                                             PacketLengthType, ...
                                                                                             PacketLength, ...
                                                                                             tm_mode.N_eff_TX, ...
                                                                                             u);
        
            % 5.2.2
            [physical_resource_mapping_STF_cell] = lib_5_physical_layer_transmissions.STF(numerology, ...
                                                                                          k_b_OCC, ...
                                                                                          tm_mode.N_eff_TX, ...
                                                                                          b);
        
            % 5.2.3
            [physical_resource_mapping_DRS_cell] = lib_5_physical_layer_transmissions.DRS(numerology, ...
                                                                                          k_b_OCC, ...
                                                                                          tm_mode.N_TS, ...
                                                                                          tm_mode.N_eff_TX, ...
                                                                                          N_PACKET_symb, ...
                                                                                          b);
        
            % 5.2.4
            [physical_resource_mapping_PCC_cell] = lib_5_physical_layer_transmissions.PCC(numerology, ...
                                                                                          k_b_OCC, ...
                                                                                          tm_mode.N_TS, ...
                                                                                          N_PACKET_symb, ...
                                                                                          physical_resource_mapping_STF_cell, ...
                                                                                          physical_resource_mapping_DRS_cell);
        
            % 5.2.5
            [physical_resource_mapping_PDC_cell, N_PDC_subc] = lib_5_physical_layer_transmissions.PDC(u, ...
                                                                                                      numerology, ...
                                                                                                      k_b_OCC, ...
                                                                                                      N_PACKET_symb, ...
                                                                                                      tm_mode.N_TS, ...
                                                                                                      tm_mode.N_eff_TX, ...
                                                                                                      physical_resource_mapping_STF_cell, ...
                                                                                                      physical_resource_mapping_DRS_cell, ...
                                                                                                      physical_resource_mapping_PCC_cell);
        
            % 5.3
            N_TB_bits = lib_5_physical_layer_transmissions.Transport_block_size(tm_mode, mcs, N_PDC_subc, Z);
        
            % save data in structure
            obj.phy_4_5.tm_mode                             = tm_mode;
            obj.phy_4_5.mcs                                 = mcs;
            obj.phy_4_5.numerology                          = numerology;
            obj.phy_4_5.T_frame                             = T_frame;
            obj.phy_4_5.N_FRAME_slot                        = N_FRAME_slot;
            obj.phy_4_5.T_slot                              = T_slot;
            obj.phy_4_5.k_b_OCC                             = k_b_OCC;
            obj.phy_4_5.N_PACKET_symb                       = N_PACKET_symb;
            obj.phy_4_5.physical_resource_mapping_STF_cell  = physical_resource_mapping_STF_cell;
            obj.phy_4_5.physical_resource_mapping_DRS_cell  = physical_resource_mapping_DRS_cell;
            obj.phy_4_5.physical_resource_mapping_PCC_cell  = physical_resource_mapping_PCC_cell;
            obj.phy_4_5.physical_resource_mapping_PDC_cell  = physical_resource_mapping_PDC_cell;
            obj.phy_4_5.N_PDC_subc                          = N_PDC_subc;
            obj.phy_4_5.N_TB_bits                           = N_TB_bits;
        
            % custom values, all starting with n_
        
            % how many gross bits can we transmit?
            obj.phy_4_5.G = tm_mode.N_SS*N_PDC_subc*mcs.N_bps;
        
            % what percentage of the spectrum do we occupy?
            obj.phy_4_5.n_spectrum_occupied = numel(k_b_OCC)/numerology.N_b_DFT;
        
            % how long is a symbol (CP included) in samples?
            obj.phy_4_5.n_T_u_symb_samples = (obj.phy_4_5.numerology.N_b_DFT*9)/8;
        
            % how long is a packet in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
            obj.phy_4_5.n_packet_samples = obj.phy_4_5.N_PACKET_symb*obj.phy_4_5.n_T_u_symb_samples;
        
            % How long are STF, DF and GI in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
            % How often does the pattern in STF repeat? (Figures 5.1-1, 5.1-2, 5.1-3)
            switch u
                case 1
                    obj.phy_4_5.n_STF_samples = (obj.phy_4_5.n_T_u_symb_samples*14)/9;
                    obj.phy_4_5.n_DF_samples = (obj.phy_4_5.N_PACKET_symb-2)*obj.phy_4_5.n_T_u_symb_samples;
                    obj.phy_4_5.n_GI_samples = (obj.phy_4_5.n_T_u_symb_samples*4)/9;
                    obj.phy_4_5.n_STF_pattern = 7;
                case {2,4}
                    obj.phy_4_5.n_STF_samples = obj.phy_4_5.n_T_u_symb_samples*2;
                    obj.phy_4_5.n_DF_samples = (obj.phy_4_5.N_PACKET_symb-3)*obj.phy_4_5.n_T_u_symb_samples;
                    obj.phy_4_5.n_GI_samples = obj.phy_4_5.n_T_u_symb_samples;
                    obj.phy_4_5.n_STF_pattern = 9;
                case 8
                    obj.phy_4_5.n_STF_samples = obj.phy_4_5.n_T_u_symb_samples*2;
                    obj.phy_4_5.n_DF_samples = (obj.phy_4_5.N_PACKET_symb-4)*obj.phy_4_5.n_T_u_symb_samples;
                    obj.phy_4_5.n_GI_samples = obj.phy_4_5.n_T_u_symb_samples*2;
                    obj.phy_4_5.n_STF_pattern = 9;                  
            end
        
            assert(obj.phy_4_5.n_packet_samples == obj.phy_4_5.n_STF_samples + obj.phy_4_5.n_DF_samples + obj.phy_4_5.n_GI_samples);
        
            obj.phy_4_5.n_pcc_bits_uncoded = 196;
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

            mcs                 = obj.phy_4_5.mcs;

            N_b_DFT             = obj.phy_4_5.numerology.N_b_DFT;
            N_b_CP              = obj.phy_4_5.numerology.N_b_CP;

            N_PACKET_symb       = obj.phy_4_5.N_PACKET_symb;
            k_b_OCC             = obj.phy_4_5.k_b_OCC;

            G                   = obj.phy_4_5.G;

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

            %% clause 7, based on the generic procedures of clause 6

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
                                                                            G, ...
                                                                            Z, ...
                                                                            network_id, ...
                                                                            PLCF_type, ...
                                                                            rv, ...
                                                                            mcs, ...
                                                                            N_SS);
                                                                        
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
            obj.packet_data.pcc_enc_dbg = pcc_enc_dbg;
            obj.packet_data.x_PDC = x_PDC;
            obj.packet_data.pdc_enc_dbg = pdc_enc_dbg;
            obj.packet_data.x_PCC_ss = x_PCC_ss;
            obj.packet_data.x_PDC_ss = x_PDC_ss;
            obj.packet_data.y_PCC_ts = y_PCC_ts;
            obj.packet_data.y_PDC_ts = y_PDC_ts;
            obj.packet_data.antenna_streams_mapped = antenna_streams_mapped;
            obj.packet_data.samples_antenna_tx = samples_antenna_tx;

            if verbosity_ > 0
                obj.plot_resource_allocation_in_frequency_domain();
                obj.plot_resource_allocation_in_time_domain();
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

        function [] = plot_resource_allocation_in_frequency_domain(obj)
            % create matrix with all components
            [~, mat_STF_DRS_PCC_PDC_all_streams] = lib_util.matrix_STF_DRS_PCC_PDC(obj.phy_4_5.numerology.N_b_DFT, ...
                                                                                   obj.phy_4_5.N_PACKET_symb, ...
                                                                                   obj.phy_4_5.tm_mode.N_TS, ...
                                                                                   obj.phy_4_5.tm_mode.N_SS, ...
                                                                                   obj.phy_4_5.physical_resource_mapping_STF_cell, ...
                                                                                   obj.phy_4_5.physical_resource_mapping_DRS_cell, ...
                                                                                   obj.phy_4_5.physical_resource_mapping_PCC_cell, ...
                                                                                   obj.phy_4_5.physical_resource_mapping_PDC_cell);
        
            figure()
            clf()
            imagesc(mat_STF_DRS_PCC_PDC_all_streams);
            title('TX Resource Allocation in Frequency Domain');
            set(gca, 'YTick', 0:2:obj.phy_4_5.numerology.N_b_DFT, 'YTickLabel', obj.phy_4_5.numerology.N_b_DFT/2-(0:2:obj.phy_4_5.numerology.N_b_DFT))
            ylabel('Subcarrier Index');
            xlabel('OFDM symbol index');
            axis image
            axis ij
            colormap jet
            clim([0 4])
        end

        % This method must only be called after a packet was generated, otherwise the time domain signal is missing.
        function [] = plot_resource_allocation_in_time_domain(obj)
            figure()
            clf()

            for i=1:1:size(obj.packet_data.samples_antenna_tx, 2)
                subplot(size(obj.packet_data.samples_antenna_tx, 2), 1, i)
                plot(abs(obj.packet_data.samples_antenna_tx(:,i)));
                hold on
                yline(rms(obj.packet_data.samples_antenna_tx(:,i)), "LineWidth", 2);
                legend("IQ", "RMS");
                title("TX Resource Allocation in Time Domain of Antenna " + num2str(i));
            end

            ylabel('Absolute Amplitude');
            xlabel('Samples Index');
            ylim([0 3])
        end
    end
end
