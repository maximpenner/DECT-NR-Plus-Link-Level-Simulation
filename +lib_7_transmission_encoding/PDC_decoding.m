function [tb_bits_recovered, PDC_HARQ_buf_report, pdc_dec_dbg] = PDC_decoding(x_PDC, ...
                                                                              N_TB_bits, ...
                                                                              Z, ...
                                                                              network_id, ...
                                                                              PLCF_type, ...
                                                                              rv, ...
                                                                              modulation, ...
                                                                              HARQ_buf)
    %%
    d = lteSymbolDemodulate(x_PDC, modulation, 'Soft');
    
    % required for comparing BER
    d_hard = lteSymbolDemodulate(x_PDC, modulation, 'Hard');
    
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
    [seq,~] = ltePRBS(g_init, numel(d),'signed');
    f = d.*seq;    

    %% rate matching
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

    cbsbuffers = HARQ_buf;

    d_turbo = lteRateRecoverTurbo(f, N_TB_bits, rv, chs, cbsbuffers);
    
    PDC_HARQ_buf_report = d_turbo; % the output of lteRateRecoverTurbo() can be passed on as a cbsbuffers in the next call, will be cleared if CRC is correct

    %%
    NTurboDecIts = 5;
    c = lteTurboDecode(d_turbo,NTurboDecIts);

    %%
    % depends on Z
    if Z == 2048
        [b,segErr] = lib_6_generic_procedures.lib_cc_rm_i.Code_block_desegmentation_Z_2048(c, (N_TB_bits+24));
    elseif Z == 6144
        [b,segErr] = lteCodeBlockDesegment(c,(N_TB_bits+24));
    else
        error('Z must be 2048 or 6144, it is %d.',Z);
    end

    %%
    [a,crcError] = lteCRCDecode(b,'24A');
    
    %% crcError is a logical values that indicates if the crc across the data and the crc appended are the same
    if crcError == 0
        tb_bits_recovered = a;
        PDC_HARQ_buf_report = [];   % bits are correct, no need fo buffer in next call
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
