function [tb_bits_recovered, PDC_HARQ_buf_report, pdc_dec_dbg] = PDC_decoding(x_PDC, ...
                                                                              N_TB_bits, ...
                                                                              Z, ...
                                                                              network_id, ...
                                                                              PLCF_type, ...
                                                                              rv, ...
                                                                              mcs, ...
                                                                              N_SS, ...
                                                                              HARQ_buf)
    %% 7.6.7 Symbol mapping
    d = lteSymbolDemodulate(x_PDC, mcs.modulation0, 'Soft');
    
    % required for comparing BER
    d_hard = lteSymbolDemodulate(x_PDC, mcs.modulation0, 'Hard');
    
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
    [seq,~] = ltePRBS(g_init, numel(d),'signed');
    f = d.*seq;    

    %% 7.6.4 Channel coding & rate matching and 7.6.5 Code block concatenation
    chs.Modulation = mcs.modulation0;
    chs.NLayers = N_SS;
    %chs.TxScheme
    %chs.NIR
    %chs.NSoftbits
    %chs.DuplexMode
    %chs.TDDConfig

    cbsbuffers = HARQ_buf;

    if Z == 2048
        G = numel(f);
        N_L = N_SS;
        Q_m = mcs.N_bps;
        d_turbo = dectnrp_6_generic_procedures.lib_cc_rm_i.lteRateRecoverTurbo(f, N_TB_bits, rv, chs, cbsbuffers, G, N_L, Q_m);
    elseif Z == 6144
        d_turbo = lteRateRecoverTurbo(f, N_TB_bits, rv, chs, cbsbuffers);
    else
        error('Z must be 2048 or 6144, it is %d.',Z);
    end
    
    % the output of lteRateRecoverTurbo() can be passed on as a cbsbuffers in the next call, will be cleared if CRC is correct
    PDC_HARQ_buf_report = d_turbo;

    NTurboDecIts = 5;
    c = lteTurboDecode(d_turbo, NTurboDecIts);

    %% 7.6.3 Code block segmentation
    if Z == 2048
        [b,segErr] = dectnrp_6_generic_procedures.lib_cc_rm_i.code_block_desegmentation(c, N_TB_bits+24);
    elseif Z == 6144
        [b,segErr] = lteCodeBlockDesegment(c, N_TB_bits+24);
    else
        error('Z must be 2048 or 6144, it is %d.',Z);
    end

    %% 7.6.2 CRC calculation
    [a,crcError] = lteCRCDecode(b,'24A');
    
    if crcError == 0
        tb_bits_recovered = a;

        % bits are correct, no need for buffer in next call
        PDC_HARQ_buf_report = [];
    else
        tb_bits_recovered = [];
    end    

    %% create an output structure for debugging
    pdc_dec_dbg.a = a;
    pdc_dec_dbg.b = b;
    pdc_dec_dbg.c = c;
    pdc_dec_dbg.d_turbo = d_turbo;
    pdc_dec_dbg.f = f;
    pdc_dec_dbg.d = d;
    pdc_dec_dbg.d_hard = d_hard;
    pdc_dec_dbg.n_code_block_error = sum(segErr~=0);
end
