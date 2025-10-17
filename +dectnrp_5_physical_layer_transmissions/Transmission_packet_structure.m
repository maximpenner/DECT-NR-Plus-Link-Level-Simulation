function [N_PACKET_symb] = Transmission_packet_structure(numerology, PacketLengthType, PacketLength, N_eff_TX, u)

    switch PacketLengthType
        case 0
            N_PACKET_symb = PacketLength * numerology.N_SLOT_u_symb / numerology.N_SLOT_u_subslot;
        case 1
            N_PACKET_symb = PacketLength * numerology.N_SLOT_u_symb;
        otherwise
            error('PacketLengthType is %d, but it must be 0 or 1.', PacketLengthType);
    end

    % minimum length to accommodate two DRS symbols per transmit stream
    assert(~(N_eff_TX == 4 && N_PACKET_symb < 15), ...
           'For N_eff_TX == 4 at least 15 symbols are required, see 5.1 Transmission packet structure.');

    % collision between first zero-symbol in GI and final DRS symbol
    if u == 8 && N_eff_TX == 8
        if N_PACKET_symb < 20 || mod(N_PACKET_symb, 10) ~= 0
            error('For u=8 and N_eff_TX=8 at least 20 symbols are required and N_PACKET_symb must be a multiple of 10, see 5.1 Transmission packet structure.');
        end
    end

    % Figures 4.5-2 and 4.5-3
    tx_packet_len = N_PACKET_symb*numerology.T_u_symb;
    tx_packet_len_STF = 14/9*numerology.T_u_symb;
    tx_packet_len_DF = (N_PACKET_symb-2)*numerology.T_u_symb;
    tx_packet_len_GI = 4/9*numerology.T_u_symb;

    assert(round(tx_packet_len) == round(tx_packet_len_STF+tx_packet_len_DF+tx_packet_len_GI));
end
