function [plcf_bits_recovered, ...
          PCC_HARQ_buf_40_report, ...
          PCC_HARQ_buf_80_report, ...
          CL_report, ...
          BF_report, ...
          pcc_dec_dbg] = PCC_decoding(x_PCC, PLCF_type, HARQ_buf_40, HARQ_buf_80)

    % notation following Figure 7.5.1-1: Physical control channel encoding

    %% 7.5.5 Symbol mapping
    % d has a positive sign for bit 1 and a negative sign for bit 0
    d = lteSymbolDemodulate(x_PCC, 'QPSK', 'Soft');
    d_hard = lteSymbolDemodulate(x_PCC, 'QPSK', 'Hard');

    %% 7.5.4 Scrambling
    % we must use soft decoding
    [seq,~] = ltePRBS(hex2dec('0x44454354'), 98*2, 'signed');
    
    % if seq is logical 1, switch the sign, if seq is logical 0, keep the sign
    e = d.*seq;   

    %% blindly determine PLCF size
    
    % According to 7.5.1, we blind decode both possible lengths and check which CRC is correct.
    % Note that PCC decoding is often ambiguous, i.e. we can have a correct CRC for 40 bit even if 80 were transmitted.
    % In a real receiver, we could easily avoid that by also checking the PCC content.
    % In this simulation, the PLCF content is random. So instead we assume we know the exact type.

    assert(PLCF_type == 1 || PLCF_type == 2);

    n_bits = 40;

    if PLCF_type == 2
        n_bits = 80;
    end

    % provide default values
    PCC_HARQ_buf_40_report = [];
    PCC_HARQ_buf_80_report = [];

    %% 7.5.3 Channel coding & rate matching
    chs.Modulation = 'QPSK';
    %chs.NLayers
    %chs.TxScheme
    %chs.NIR
    %chs.NSoftbits
    %chs.DuplexMode
    %chs.TDDConfig

    if n_bits == 40
        cbsbuffers = HARQ_buf_40;
    else
        cbsbuffers = HARQ_buf_80;
    end

    % Why (n_bits-8)?
    %
    % PCC user bits is either n_bits=40 or n_bits=80 bits plus 16 bit CRC.
    % lteRateRecoverTurbo, however, assumes a 24 bit CRC.
    % So we pretend the last 8 bits of the PCC user bits belong to the CRC.
    % This way it looks like we passed on either 32 or 72 bit.
    rv = 0;
    d_turbo = lteRateRecoverTurbo(e, (n_bits-8), rv, chs, cbsbuffers);
    
    if n_bits == 40
        PCC_HARQ_buf_40_report = d_turbo;
    else
        PCC_HARQ_buf_80_report = d_turbo;
    end        

    NTurboDecIts = 5;
    c = lteTurboDecode(d_turbo,NTurboDecIts);

    %% 7.5.2 CRC calculation
    mask_CL = logical(repmat([0, 1], 1, 8));
    mask_BF = logical(repmat([1, 0], 1, 8));
    mask_CL_BF = xor(mask_CL, mask_BF);
    
    [a,crcError] = lteCRCDecode(cell2mat(c),'16');
    
    if crcError == 0
        plcf_bits_recovered = a;
        CL_report = false;
        BF_report = false;
    elseif xor(de2bi(crcError,16,'left-msb'), mask_CL_BF) == 0
        plcf_bits_recovered = a;
        CL_report = true;
        BF_report = true;
    elseif xor(de2bi(crcError,16,'left-msb'), mask_BF) == 0
        plcf_bits_recovered = a;
        CL_report = false;
        BF_report = true;
    elseif xor(de2bi(crcError,16,'left-msb'), mask_CL) == 0
        plcf_bits_recovered = a;
        CL_report = true;
        BF_report = false;
    else
        plcf_bits_recovered = [];
        CL_report = [];
        BF_report = [];
    end
    
    %% create an output structure for debugging
    pcc_dec_dbg.d = d;
    pcc_dec_dbg.d_hard = d_hard;
    pcc_dec_dbg.e = e;
    pcc_dec_dbg.d_turbo = d_turbo;
    pcc_dec_dbg.c = c;
    pcc_dec_dbg.a = a;
end
