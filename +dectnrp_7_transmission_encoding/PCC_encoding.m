function [x_PCC, pcc_enc_dbg] = PCC_encoding(plcf_bits, CL, is_precoding_identity_matrix)

    % notation following Figure 7.5.1-1: Physical control channel encoding

    a = plcf_bits;

    %% 7.5.2 CRC calculation
    mask_CL = logical(repmat([0, 1], 1, 8));
    mask_BF = logical(repmat([1, 0], 1, 8));
    mask_CL_BF = xor(mask_CL, mask_BF);
    
    if CL == true && is_precoding_identity_matrix == true
        mask = mask_CL;
    elseif CL == false && is_precoding_identity_matrix == false
        mask = mask_BF;
    elseif CL == true && is_precoding_identity_matrix == false
        mask = mask_CL_BF;
    else
        mask = false(1,16);
    end
    
    % mask is a logical with the leftest value being the msb, mask input to lteCRCEncode() must be an integer
    c = lteCRCEncode(a,'16', bi2de(mask, 'left-msb'));
    
    %% 7.5.3 Channel coding & rate matching
    
    % code block segmentation not applied to PCC as it is too short
    d_turbo = lteTurboEncode(c);
    
    % rate matching structure
    chs.Modulation = 'QPSK';
    %chs.NLayers
    %chs.TxScheme
    %chs.NIR
    %chs.NSoftbits
    %chs.DuplexMode
    %chs.TDDConfig

    % rv is always 0
    rv = 0;
    e = lteRateMatchTurbo(d_turbo, 98*2, rv, chs);
    
    %% 7.5.4 Scrambling
    [seq,~] = ltePRBS(hex2dec('0x44454354'), 98*2);
    d = mod(e + int8(seq), 2);

    %% 7.5.5 Symbol mapping
    x_PCC = lteSymbolModulate(d, 'QPSK');

    %% create an output structure for debugging
    pcc_enc_dbg.a = a;
    pcc_enc_dbg.c = c;
    pcc_enc_dbg.d_turbo = d_turbo;
    pcc_enc_dbg.e = e;
    pcc_enc_dbg.d = d;
end
