function [phy_4_5] = run_chapter_4_5(verbose, mac_meta)

    % This function is called in the constructor. It basically does all the calculations of chapter 4 and 5.
    % All data is saved in the structure phy_4_5.

    %% for the purpose of readability, extract all variables that are required from mac_meta
    u                   = mac_meta.u;
    b                   = mac_meta.b;
    PacketLengthType    = mac_meta.PacketLengthType;
    PacketLength        = mac_meta.PacketLength;
    tm_mode_0_to_11     = mac_meta.tm_mode_0_to_11;
    mcs_index           = mac_meta.mcs_index;
    Z                   = mac_meta.Z;

    %% start generating the frame structure with the functions provided in the technical specification (TS)

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
    N_PACKET_symb = lib_5_physical_layer_transmissions.Transmission_packet_structure(numerology, PacketLengthType, PacketLength, tm_mode.N_eff_TX, u);

    % 5.2.2
    [physical_resource_mapping_STF_cell] = lib_5_physical_layer_transmissions.STF(numerology, k_b_OCC, tm_mode.N_eff_TX, b);

    % 5.2.3
    [physical_resource_mapping_DRS_cell] = lib_5_physical_layer_transmissions.DRS(numerology, k_b_OCC, tm_mode.N_TS, tm_mode.N_eff_TX, N_PACKET_symb, b);

    % 5.2.4
    [physical_resource_mapping_PCC_cell] = lib_5_physical_layer_transmissions.PCC(numerology, k_b_OCC, tm_mode.N_TS, N_PACKET_symb,...
                                                                                    physical_resource_mapping_STF_cell,...
                                                                                    physical_resource_mapping_DRS_cell);

    % 5.2.5
    [physical_resource_mapping_PDC_cell, N_PDC_subc] = lib_5_physical_layer_transmissions.PDC(u, numerology, k_b_OCC, N_PACKET_symb, tm_mode.N_TS, tm_mode.N_eff_TX,...
                                                                                                physical_resource_mapping_STF_cell,...
                                                                                                physical_resource_mapping_DRS_cell,...
                                                                                                physical_resource_mapping_PCC_cell);

    % 5.3
    N_TB_bits = lib_5_physical_layer_transmissions.Transport_block_size(tm_mode, mcs, N_PDC_subc, Z);

    % save data in structure
    phy_4_5.tm_mode                             = tm_mode;
    phy_4_5.mcs                                 = mcs;
    phy_4_5.numerology                          = numerology;
    phy_4_5.T_frame                             = T_frame;
    phy_4_5.N_FRAME_slot                        = N_FRAME_slot;
    phy_4_5.T_slot                              = T_slot;
    phy_4_5.k_b_OCC                             = k_b_OCC;
    phy_4_5.N_PACKET_symb                       = N_PACKET_symb;
    phy_4_5.physical_resource_mapping_STF_cell  = physical_resource_mapping_STF_cell;
    phy_4_5.physical_resource_mapping_DRS_cell  = physical_resource_mapping_DRS_cell;
    phy_4_5.physical_resource_mapping_PCC_cell  = physical_resource_mapping_PCC_cell;
    phy_4_5.physical_resource_mapping_PDC_cell  = physical_resource_mapping_PDC_cell;
    phy_4_5.N_PDC_subc                          = N_PDC_subc;
    phy_4_5.N_TB_bits                           = N_TB_bits;

    % custom values, all starting with n_

    % how many gross bits can we transmit?
    phy_4_5.n_total_bits                        = tm_mode.N_SS*N_PDC_subc*mcs.N_bps;

    % what percentage of the spectrum do we occupy?
    phy_4_5.n_spectrum_occupied                 = numel(k_b_OCC)/numerology.N_b_DFT;

    % how long is a symbol (CP included) in samples?
    phy_4_5.n_T_u_symb_samples                  = (phy_4_5.numerology.N_b_DFT*9)/8;

    % how long is a packet in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
    phy_4_5.n_packet_samples                    = phy_4_5.N_PACKET_symb*phy_4_5.n_T_u_symb_samples;

    % How long are STF, DF and GI in samples? (Figures 5.1-1, 5.1-2, 5.1-3)
    % How often does the pattern in STF repeat? (Figures 5.1-1, 5.1-2, 5.1-3)
    switch u
        case 1
            phy_4_5.n_STF_samples               = (phy_4_5.n_T_u_symb_samples*14)/9;
            phy_4_5.n_DF_samples                = (phy_4_5.N_PACKET_symb-2)*phy_4_5.n_T_u_symb_samples;
            phy_4_5.n_GI_samples                = (phy_4_5.n_T_u_symb_samples*4)/9;
            phy_4_5.n_STF_pattern               = 7;
        case {2,4}
            phy_4_5.n_STF_samples               = phy_4_5.n_T_u_symb_samples*2;
            phy_4_5.n_DF_samples                = (phy_4_5.N_PACKET_symb-3)*phy_4_5.n_T_u_symb_samples;
            phy_4_5.n_GI_samples                = phy_4_5.n_T_u_symb_samples;
            phy_4_5.n_STF_pattern               = 9;
        case 8
            phy_4_5.n_STF_samples               = phy_4_5.n_T_u_symb_samples*2;
            phy_4_5.n_DF_samples                = (phy_4_5.N_PACKET_symb-4)*phy_4_5.n_T_u_symb_samples;
            phy_4_5.n_GI_samples                = phy_4_5.n_T_u_symb_samples*2;
            phy_4_5.n_STF_pattern               = 9;                  
    end

    % sanity check
    if phy_4_5.n_packet_samples ~= phy_4_5.n_STF_samples + phy_4_5.n_DF_samples + phy_4_5.n_GI_samples
        error('Length of packet is not the sum of STF, DF and GI.');
    end

    % debugging
    if verbose > 1
        lib_util.plot_resource_mapping(numerology, N_PACKET_symb, tm_mode.N_TS, tm_mode.N_SS,...
                                        physical_resource_mapping_STF_cell,...
                                        physical_resource_mapping_DRS_cell,...
                                        physical_resource_mapping_PCC_cell,...
                                        physical_resource_mapping_PDC_cell);
    end
end

