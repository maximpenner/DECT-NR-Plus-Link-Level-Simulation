function [PCC_user_bits,...
            PCC_harq_buf_40_report,...
            PCC_harq_buf_80_report,...
            CL_report,...
            BF_report,...
            PCC_report,...
            pcc_dec_dbg] = PCC_decoding(x_PCC, harq_buf_40, harq_buf_80)

    %%
    % d has a positive sign for bit 1 and a negative sign for bit 0
    d = lteSymbolDemodulate(x_PCC, 'QPSK', 'Soft');
    d_hard = lteSymbolDemodulate(x_PCC, 'QPSK', 'Hard');

    %% descramble
    % we must use softdecoding
    [seq,~] = ltePRBS(hex2dec('0x44454354'), 98*2, 'signed');
    
    % if seq is logical 1, switch the sign, if seq is logical 0, keep the sign
    e = d.*seq;   

    %% according to 7.5.1 we blind decode both possible lengths and check which crc is correct
    for n_bits = [40, 80]
        %% channel coding and rate matching
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

        if n_bits == 40
            cbsbuffers = harq_buf_40;
        else
            cbsbuffers = harq_buf_80;
        end
        %
        % Why (n_bits-8)?
        %
        % PCC user bits is either n_bits=40 or n_bits=80 bits plus 16 bit CRC.
        % lteRateRecoverTurbo, however, assumes a 24 bit CRC.
        % So we pretend the last 8 bits of the PCC user bits belong to the CRC.
        % This way it looks like we passed on either 32 or 72 bit.
        %
        rv = 0;     % see 7.5.3
        d_turbo = lteRateRecoverTurbo(e, (n_bits-8), rv, chs, cbsbuffers);
        
        if n_bits == 40
            PCC_harq_buf_40_report = d_turbo;
        else
            PCC_harq_buf_80_report = d_turbo;
        end        

        % channel decoding
        NTurboDecIts = 5;
        c = lteTurboDecode(d_turbo,NTurboDecIts);

        %% code block segmentation not applied to PCC as it is to short
        %

        %% crc calculation
        % hexToBinaryVector requires system control toolbox
        %mask_CL = hexToBinaryVector('0x5555',16);   % 0101010101010101 = (21845)_10, bi2de(mask_CL,'left-msb')=21845
        %mask_BF = hexToBinaryVector('0xAAAA',16);   % 1010101010101010 = (43690)_10
        mask_CL = logical(repmat([0, 1], 1, 8));
        mask_BF = logical(repmat([1, 0], 1, 8));
        mask_CL_BF = xor(mask_CL, mask_BF);
        
        [a,crcError] = lteCRCDecode(cell2mat(c),'16');
        
        if crcError == 0
            PCC_user_bits = a;
            CL_report = false;
            BF_report = false;
        elseif xor(de2bi(crcError,16,'left-msb'), mask_CL_BF) == 0
            PCC_user_bits = a;
            CL_report = true;
            BF_report = true;
        elseif xor(de2bi(crcError,16,'left-msb'), mask_BF) == 0
            PCC_user_bits = a;
            CL_report = false;
            BF_report = true;
        elseif xor(de2bi(crcError,16,'left-msb'), mask_CL) == 0
            PCC_user_bits = a;
            CL_report = true;
            BF_report = false;
        else
            PCC_user_bits = [];
            CL_report = [];
            BF_report = [];
            PCC_report = [];
        end
        
        %% create an output structure for debugging
        pcc_dec_dbg.d = d;
        pcc_dec_dbg.d_hard = d_hard;
        pcc_dec_dbg.e = e;
        pcc_dec_dbg.d_turbo = d_turbo;
        pcc_dec_dbg.c = c;
        pcc_dec_dbg.a = a;
        
        %% leave if any crc for 40 bit was correct
        if numel(PCC_user_bits) > 0
            
            PCC_harq_buf_40_report = [];
            PCC_harq_buf_80_report = [];
            
            % TODO: extract all relevant information in PCC
            PCC_report = [];
            
            break;
        end
    end
end

