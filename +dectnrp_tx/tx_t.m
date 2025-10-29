classdef tx_t < matlab.mixin.Copyable
    
    properties
        config;

        % The following structure contains variables that are derived from the above minimal set of variables in config.
        % They are described in clause 4 and 5 of ETSI TS 103 636-3.
        derived;

        % intermediate results during packet encoding
        packet_data;
    end
    
    methods
        function obj = tx_t(config)
            assert(isa(config, "dectnrp_tx.config_t"));
            
            obj.config = config;
            obj.set_derived();
            obj.packet_data = [];
        end
        
        function [samples_antenna_tx] = generate_packet(obj, plcf_bits, tb_bits)
            
            %% for the purpose of readability, read all variables that are necessary at this stage

            u               = obj.config.u;
            b               = obj.config.b;
            Z               = obj.config.Z;
            codebook_index  = obj.config.codebook_index;
            network_id      = obj.config.network_id;
            PLCF_type       = obj.config.PLCF_type;
            rv              = obj.config.rv;
            oversampling    = obj.config.oversampling;
            verbosity       = obj.config.verbosity;

            mode_0_to_11    = obj.derived.tm_mode.mode_0_to_11;
            N_SS            = obj.derived.tm_mode.N_SS;
            CL              = obj.derived.tm_mode.CL;
            N_TS            = obj.derived.tm_mode.N_TS;
            N_TX            = obj.derived.tm_mode.N_TX;
            N_eff_TX        = obj.derived.tm_mode.N_eff_TX;

            mcs             = obj.derived.mcs;

            N_b_DFT         = obj.derived.numerology.N_b_DFT;
            N_b_CP          = obj.derived.numerology.N_b_CP;

            N_PACKET_symb   = obj.derived.N_PACKET_symb;
            k_b_OCC         = obj.derived.k_b_OCC;

            G               = obj.derived.G;

            physical_resource_mapping_PCC_cell = obj.derived.physical_resource_mapping_PCC_cell;
            physical_resource_mapping_PDC_cell = obj.derived.physical_resource_mapping_PDC_cell;
            physical_resource_mapping_STF_cell = obj.derived.physical_resource_mapping_STF_cell;
            physical_resource_mapping_DRS_cell = obj.derived.physical_resource_mapping_DRS_cell;

            %% clause 7, based on the generic procedures of clause 6

            % The receiver needs to know if signal is beamformed or not for channel sounding purposes.
            % 7.2
            if mode_0_to_11 == 0
                is_precoding_identity_matrix = true;
            else
                if codebook_index == 0
                    is_precoding_identity_matrix = true;
                else
                    is_precoding_identity_matrix = false;
                end
            end
            
            % First, we calculate complex samples for PCC and PDC.
            % pcc_enc_dbg and pdc_enc_dbg contain debugging information.
            [x_PCC, pcc_enc_dbg] = dectnrp_7_transmission_encoding.PCC_encoding(plcf_bits, ...
                                                                                CL, ...
                                                                                is_precoding_identity_matrix);
            [x_PDC, pdc_enc_dbg] = dectnrp_7_transmission_encoding.PDC_encoding(tb_bits, ...
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
                x_PDC_ss = dectnrp_6_generic_procedures.Spatial_Multiplexing(x_PDC, N_SS);
            else
                x_PDC_ss = {x_PDC};
            end
                                                                        
            % Transmit diversity precoding, switching to transmit streams (ts).
            % Always applied to PCC if more than one transmit stream is used.
            % Applied to PDC only if we are in mode 1, 5 or 10 according to Table 7.2-1 with N_SS = 1.
            % Otherwise, we just directly switch from spatial streams to transmit streams according to 6.3.3.1.
            if N_TS > 1
                y_PCC_ts = dectnrp_6_generic_procedures.Transmit_diversity_precoding(x_PCC_ss, N_TS);
            else
                y_PCC_ts = x_PCC_ss;
            end
            if ismember(mode_0_to_11, [1,5,10]) == true
                y_PDC_ts = dectnrp_6_generic_procedures.Transmit_diversity_precoding(x_PDC_ss, N_TS);
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
            transmit_streams = dectnrp_7_transmission_encoding.subcarrier_mapping_STF(transmit_streams, physical_resource_mapping_STF_cell);
            transmit_streams = dectnrp_7_transmission_encoding.subcarrier_mapping_DRS(transmit_streams, physical_resource_mapping_DRS_cell);            
            transmit_streams = dectnrp_7_transmission_encoding.subcarrier_mapping_PCC(transmit_streams, physical_resource_mapping_PCC_cell, y_PCC_ts);
            transmit_streams = dectnrp_7_transmission_encoding.subcarrier_mapping_PDC(transmit_streams, physical_resource_mapping_PDC_cell, y_PDC_ts);
            
            % Beamforming (N_eff_TX many), remember N_eff_TX = N_TS
            antenna_streams_mapped = dectnrp_6_generic_procedures.Beamforming(transmit_streams, N_TX, codebook_index);

            % switch to time domain
            samples_antenna_tx = dectnrp_6_generic_procedures.ofdm_signal_generation_Cyclic_prefix_insertion(antenna_streams_mapped, ...
                                                                                                             k_b_OCC, ...
                                                                                                             N_PACKET_symb, ...
                                                                                                             N_TX, ...
                                                                                                             N_eff_TX, ...
                                                                                                             N_b_DFT, ...
                                                                                                             u, ...
                                                                                                             N_b_CP, ...
                                                                                                             oversampling);

            % apply STF cover sequence
            samples_antenna_tx = dectnrp_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_tx, u, b*oversampling);

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

            if verbosity > 0
                obj.plot_resource_allocation_in_frequency_domain();
                obj.plot_resource_allocation_in_time_domain();
            end
        end

        function [samples_antenna_tx] = generate_random_packet(obj)

            % PLCF
            if obj.config.PLCF_type == 1
                plcf_bits = randi([0 1], 40, 1);
            elseif obj.config.PLCF_type == 2
                plcf_bits = randi([0 1], 80, 1);
            else
                error("undefined PLCF type");
            end
            
            % transport block
            tb_bits = randi([0 1], obj.derived.N_TB_bits, 1);

            samples_antenna_tx = obj.generate_packet(plcf_bits, tb_bits);
        end

        function [] = plot_resource_allocation_in_frequency_domain(obj)
            % create matrix with all components
            [~, mat_STF_DRS_PCC_PDC_all_streams] = dectnrp_5_physical_layer_transmissions.matrix_STF_DRS_PCC_PDC(obj.derived.numerology.N_b_DFT, ...
                                                                                                                 obj.derived.N_PACKET_symb, ...
                                                                                                                 obj.derived.tm_mode.N_TS, ...
                                                                                                                 obj.derived.tm_mode.N_SS, ...
                                                                                                                 obj.derived.physical_resource_mapping_STF_cell, ...
                                                                                                                 obj.derived.physical_resource_mapping_DRS_cell, ...
                                                                                                                 obj.derived.physical_resource_mapping_PCC_cell, ...
                                                                                                                 obj.derived.physical_resource_mapping_PDC_cell);
        
            figure()
            clf()
            imagesc(mat_STF_DRS_PCC_PDC_all_streams);
            title('TX Resource Allocation in Frequency Domain');
            set(gca, 'YTick', 0:2:obj.derived.numerology.N_b_DFT, 'YTickLabel', obj.derived.numerology.N_b_DFT/2-(0:2:obj.derived.numerology.N_b_DFT))
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

    methods (Hidden = true)
        % This method is called in the constructor. It basically does all the calculations of clause 4 and 5.
        function [] = set_derived(obj)
            % clause 7.2
            tm_mode = dectnrp_7_transmission_encoding.transmission_modes(obj.config.tm_mode_0_to_11);
        
            % Annex A
            mcs = dectnrp_Annex_A.modulation_and_coding_scheme(obj.config.mcs_index);
        
            % clause 4.3
            numerology = dectnrp_4_physical_layer_principles.numerologies(obj.config.u, obj.config.b);
        
            % clause 4.4
            [T_frame, N_FRAME_slot, T_slot] = dectnrp_4_physical_layer_principles.frame_structure();
        
            % clause 4.5
            k_b_OCC = dectnrp_4_physical_layer_principles.physical_resources(numerology.N_b_OCC);
        
            % clause 5.1
            N_PACKET_symb = dectnrp_5_physical_layer_transmissions.Transmission_packet_structure(numerology, ...
                                                                                                 obj.config.PacketLengthType, ...
                                                                                                 obj.config.PacketLength, ...
                                                                                                 tm_mode.N_eff_TX, ...
                                                                                                 obj.config.u);
            
            % clause 5.2.2
            [physical_resource_mapping_STF_cell] = dectnrp_5_physical_layer_transmissions.STF(numerology, ...
                                                                                              k_b_OCC, ...
                                                                                              tm_mode.N_eff_TX, ...
                                                                                              obj.config.b);
        
            % clause 5.2.3
            [physical_resource_mapping_DRS_cell] = dectnrp_5_physical_layer_transmissions.DRS(numerology, ...
                                                                                              k_b_OCC, ...
                                                                                              tm_mode.N_TS, ...
                                                                                              tm_mode.N_eff_TX, ...
                                                                                              N_PACKET_symb, ...
                                                                                              obj.config.b);
        
            % clause 5.2.4
            [physical_resource_mapping_PCC_cell] = dectnrp_5_physical_layer_transmissions.PCC(numerology, ...
                                                                                              k_b_OCC, ...
                                                                                              tm_mode.N_TS, ...
                                                                                              N_PACKET_symb, ...
                                                                                              physical_resource_mapping_STF_cell, ...
                                                                                              physical_resource_mapping_DRS_cell);
        
            % clause 5.2.5
            [physical_resource_mapping_PDC_cell, N_PDC_subc] = dectnrp_5_physical_layer_transmissions.PDC(obj.config.u, ...
                                                                                                          numerology, ...
                                                                                                          k_b_OCC, ...
                                                                                                          N_PACKET_symb, ...
                                                                                                          tm_mode.N_TS, ...
                                                                                                          tm_mode.N_eff_TX, ...
                                                                                                          physical_resource_mapping_STF_cell, ...
                                                                                                          physical_resource_mapping_DRS_cell, ...
                                                                                                          physical_resource_mapping_PCC_cell);
        
            % clause 5.3
            N_TB_bits = dectnrp_5_physical_layer_transmissions.Transport_block_size(tm_mode, mcs, N_PDC_subc, obj.config.Z);

            % clause 7.6.6
            G = tm_mode.N_SS*N_PDC_subc*mcs.N_bps;
        
            % save data in structure
            obj.derived.tm_mode                             = tm_mode;
            obj.derived.mcs                                 = mcs;
            obj.derived.numerology                          = numerology;
            obj.derived.T_frame                             = T_frame;
            obj.derived.N_FRAME_slot                        = N_FRAME_slot;
            obj.derived.T_slot                              = T_slot;
            obj.derived.k_b_OCC                             = k_b_OCC;
            obj.derived.N_PACKET_symb                       = N_PACKET_symb;
            obj.derived.physical_resource_mapping_STF_cell  = physical_resource_mapping_STF_cell;
            obj.derived.physical_resource_mapping_DRS_cell  = physical_resource_mapping_DRS_cell;
            obj.derived.physical_resource_mapping_PCC_cell  = physical_resource_mapping_PCC_cell;
            obj.derived.physical_resource_mapping_PDC_cell  = physical_resource_mapping_PDC_cell;
            obj.derived.N_PDC_subc                          = N_PDC_subc;
            obj.derived.N_TB_bits                           = N_TB_bits;
            obj.derived.G                                   = G;
        
            % custom values, all starting with n_
        
            % what percentage of the spectrum do we occupy?
            obj.derived.n_spectrum_occupied = numel(k_b_OCC)/numerology.N_b_DFT;
        
            % how long is a symbol (CP included) in samples?
            obj.derived.n_T_u_symb_samples = (obj.derived.numerology.N_b_DFT*9)/8;
        
            % how long is a packet in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
            obj.derived.n_packet_samples = obj.derived.N_PACKET_symb*obj.derived.n_T_u_symb_samples;
        
            % How long are STF, DF and GI in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
            % How often does the pattern in STF repeat? (Figures 5.1-1, 5.1-2, 5.1-3)
            switch obj.config.u
                case 1
                    obj.derived.n_STF_samples = (obj.derived.n_T_u_symb_samples*14)/9;
                    obj.derived.n_DF_samples = (obj.derived.N_PACKET_symb-2)*obj.derived.n_T_u_symb_samples;
                    obj.derived.n_GI_samples = (obj.derived.n_T_u_symb_samples*4)/9;
                    obj.derived.n_STF_pattern = 7;
                case {2,4}
                    obj.derived.n_STF_samples = obj.derived.n_T_u_symb_samples*2;
                    obj.derived.n_DF_samples = (obj.derived.N_PACKET_symb-3)*obj.derived.n_T_u_symb_samples;
                    obj.derived.n_GI_samples = obj.derived.n_T_u_symb_samples;
                    obj.derived.n_STF_pattern = 9;
                case 8
                    obj.derived.n_STF_samples = obj.derived.n_T_u_symb_samples*2;
                    obj.derived.n_DF_samples = (obj.derived.N_PACKET_symb-4)*obj.derived.n_T_u_symb_samples;
                    obj.derived.n_GI_samples = obj.derived.n_T_u_symb_samples*2;
                    obj.derived.n_STF_pattern = 9;                  
            end
        
            assert(obj.derived.n_packet_samples == obj.derived.n_STF_samples + obj.derived.n_DF_samples + obj.derived.n_GI_samples);
        
            obj.derived.n_pcc_bits_uncoded = 196;
        end
    end
end
