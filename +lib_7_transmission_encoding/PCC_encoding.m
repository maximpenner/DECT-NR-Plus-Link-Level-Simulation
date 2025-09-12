function [x_PCC, pcc_enc_dbg] = PCC_encoding(PCC_user_bits, CL, precoding_identity_matrix)
    %%
    a = PCC_user_bits;

    %% 7.5.2 crc calculation
    % convert to logical so we can xor
    % hexToBinaryVector requires system control toolbox
    %mask_CL = hexToBinaryVector('0x5555',16);   % 0101010101010101 = (21845)_10, bi2de(mask_CL,'left-msb')=21845
    %mask_BF = hexToBinaryVector('0xAAAA',16);   % 1010101010101010 = (43690)_10
    mask_CL = logical(repmat([0, 1], 1, 8));
    mask_BF = logical(repmat([1, 0], 1, 8));
    mask_CL_BF = xor(mask_CL, mask_BF);         % 1111111111111111 is this even possible?
    
    % 7.5.2.2 and 7.5.2.3
    if CL == true && precoding_identity_matrix == true
        mask = mask_CL;
    elseif CL == false && precoding_identity_matrix == false
        mask = mask_BF;
    elseif CL == true && precoding_identity_matrix == false
        mask = mask_CL_BF;
    else
        mask = false(1,16);                     % 0000000000000000
    end
    
    % mask is a logical with the leftest value being the msb, mask input to lteCRCEncode() must be an integer
    c = lteCRCEncode(a,'16', bi2de(mask,'left-msb'));
    
    %% code block segmentation not applied to PCC as it is to short
    %
    
    %% channel coding and rate matching
    d_turbo = lteTurboEncode(c);
    
    % rate matching
    chs.Modulation = 'QPSK';
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

    % rv is always 0
    rv = 0;
    e = lteRateMatchTurbo(d_turbo, 98*2, rv, chs);
    
    %% scrambling
    [seq,~] = ltePRBS(hex2dec('0x44454354'), 98*2);     % 0x44454354 = (1145389908)_10
    d = mod(e + int8(seq), 2);

    %% symbol mapping
    x_PCC = lteSymbolModulate(d,'QPSK');

    %% create an output structure for debugging
    pcc_enc_dbg.a = a;
    pcc_enc_dbg.c = c;
    pcc_enc_dbg.d_turbo = d_turbo;
    pcc_enc_dbg.e = e;
    pcc_enc_dbg.d = d;
end

