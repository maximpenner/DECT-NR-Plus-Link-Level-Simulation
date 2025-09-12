function [x_PDC, pdc_enc_dbg] = PDC_encoding(PDC_user_bits,...
                                                n_total_bits,...
                                                Z,...
                                                network_id,...
                                                PLCF_type,...
                                                rv,...
                                                modulation)
    %%
    a = PDC_user_bits;

    %%
    b = lteCRCEncode(a,'24A');

    %%
    if Z == 2048
        c = lib_6_generic_procedures.lib_cc_rm_i.Code_block_segmentation_Z_2048(b);
    elseif Z == 6144
        c = lteCodeBlockSegment(b);
    else
        error('Z must be 2048 or 6144, it is %d.',Z);
    end
    
    %%
    d_turbo = lteTurboEncode(c);
    
    %% rate matching and code block concatenation
    chs.Modulation = modulation;
    %chs.NLayers = 1;
    %chs.TxScheme = 'Port0';
    
    % determine NIR
    %NSoftbits = 25344*100;
    %M_DL_HARQ = 1;
    %M_limit = 8;
    %chs.NIR = floor(NSoftbits/min(M_DL_HARQ, M_limit));

    %chs.NSoftbits = ;
    %chs.DuplexMode = ;
    %chs.TDDConfig = ;

    f = lteRateMatchTurbo(d_turbo, n_total_bits, rv, chs);
    
    %% scrambling
    % network_id ist a 32 bit vector with network_id(1) being the MSB
    if PLCF_type == 1
        %mask_8bit = hexToBinaryVector('0x000000ff',32,'MSBFirst');
        mask_8bit = logical([zeros(1,24), ones(1,8)]);
        g_init = and(network_id,mask_8bit);
    elseif PLCF_type == 2
        network_id = [zeros(1,8) network_id(1:end-8)];
        %mask_24bit = hexToBinaryVector('0x00ffffff',32,'MSBFirst');
        mask_24bit = logical([zeros(1,8), ones(1,24)]);
        g_init = and(network_id,mask_24bit);
    end
    g_init = bi2de(g_init,'left-msb');
    [seq,~] = ltePRBS(g_init, n_total_bits);
    d = mod(f + int8(seq), 2);

    %% symbol mapping
    x_PDC = lteSymbolModulate(d,modulation);

    %% create an output structure for debugging
    pdc_enc_dbg.a = a;
    pdc_enc_dbg.b = b;
    pdc_enc_dbg.c = c;
    pdc_enc_dbg.d_turbo = d_turbo;
    pdc_enc_dbg.f = f;
    pdc_enc_dbg.d = d;
    pdc_enc_dbg.n_code_block = numel(c);
end
