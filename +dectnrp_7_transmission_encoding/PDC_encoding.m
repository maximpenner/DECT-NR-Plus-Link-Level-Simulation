function [x_PDC, pdc_enc_dbg] = PDC_encoding(tb_bits, ...
                                             G, ...
                                             Z, ...
                                             network_id, ...
                                             PLCF_type, ...
                                             rv, ...
                                             mcs, ...
                                             N_SS)
    % notation following Figure 7.6.1-1: Physical Data Channel Encoding

    a = tb_bits;

    %% 7.6.2 CRC calculation
    b = lteCRCEncode(a,'24A');

    %% 7.6.3 Code block segmentation 
    if Z == 2048
        c = dectnrp_6_generic_procedures.lib_cc_rm_i.code_block_segmentation(b);
    elseif Z == 6144
        c = lteCodeBlockSegment(b);
    else
        error('Z must be 2048 or 6144, it is %d.',Z);
    end

    assert_segmentation(numel(b), Z, c);
    
    %% 7.6.4 Channel coding & rate matching and 7.6.5 Code block concatenation
    d_turbo = lteTurboEncode(c);

    assert_turbocode(c, d_turbo);
    
    chs.Modulation = mcs.modulation0;
    chs.NLayers = N_SS;
    %chs.TxScheme
    %chs.NIR
    %chs.NSoftbits
    %chs.DuplexMode
    %chs.TDDConfig

    f = lteRateMatchTurbo(d_turbo, G, rv, chs);

    assert_ratematching(c, G, f, mcs.N_bps);
    
    %% 7.6.6 Scrambling
    % network_id is a 32 bit vector with network_id(1) being the MSB
    if PLCF_type == 1
        mask_8bit = logical([zeros(1,24), ones(1,8)]);
        g_init = and(network_id,mask_8bit);
    elseif PLCF_type == 2
        network_id = [zeros(1,8) network_id(1:end-8)];
        mask_24bit = logical([zeros(1,8), ones(1,24)]);
        g_init = and(network_id,mask_24bit);
    else
        error('PLCF_type must be 1 or 2, it is %d.', PLCF_type);
    end
    g_init = bi2de(g_init,'left-msb');
    [seq,~] = ltePRBS(g_init, G);
    d = mod(f + int8(seq), 2);

    %% 7.6.7 Symbol mapping
    x_PDC = lteSymbolModulate(d, mcs.modulation0);

    %% create an output structure for debugging
    pdc_enc_dbg.a = a;
    pdc_enc_dbg.b = b;
    pdc_enc_dbg.c = c;
    pdc_enc_dbg.d_turbo = d_turbo;
    pdc_enc_dbg.f = f;
    pdc_enc_dbg.d = d;
    pdc_enc_dbg.n_code_block = numel(c);
end

function [] = assert_segmentation(B, Z, c)

    K_r = dectnrp_6_generic_procedures.lib_cc_rm_i.get_3gpp_code_block_segment_lengths(B, Z);

    for i=1:1:numel(c)
        assert(numel(c{i}) == K_r(i));
    end
end


function [] = assert_turbocode(c, d_turbo)

    assert(numel(c) > 0);
    assert(numel(c) == numel(d_turbo));

    for i=1:1:numel(c)
        A = numel(c{i});
        B = numel(d_turbo{i});
        assert(A*3+12 == B);
    end
end

function [] = assert_ratematching(c, G, f, Q_m)

    assert(numel(c) > 0);
    assert(numel(f) == G);

    C = numel(c);
    N_L = 1;

    E_r = dectnrp_6_generic_procedures.lib_cc_rm_i.get_3gpp_encoded_code_block_segment_lengths(G, C, N_L, Q_m);

    assert(C == numel(E_r));
    assert(numel(f) == sum(E_r));
end
