function [mcs] = modulation_and_coding_scheme(MCS_index)

    switch MCS_index
        case 0
            modulation = 'BPSK';
            modulation0 = 'BPSK';
            N_bps = 1;
            R_t = 1;
            R_b = 2;
        case 1
            modulation = 'QPSK';
            modulation0 = 'QPSK';
            N_bps = 2;
            R_t = 1;
            R_b = 2;
        case 2
            modulation = 'QPSK';
            modulation0 = 'QPSK';
            N_bps = 2;
            R_t = 3;
            R_b = 4;
        case 3
            modulation = '16-QAM';
            modulation0 = '16QAM';
            N_bps = 4;
            R_t = 1;
            R_b = 2;
        case 4
            modulation = '16-QAM';
            modulation0 = '16QAM';
            N_bps = 4;
            R_t = 3;
            R_b = 4;
        case 5
            modulation = '64-QAM';
            modulation0 = '64QAM';
            N_bps = 6;
            R_t = 2;
            R_b = 3;
        case 6
            modulation = '64-QAM';
            modulation0 = '64QAM';
            N_bps = 6;
            R_t = 3;
            R_b = 4;
        case 7
            modulation = '64-QAM';
            modulation0 = '64QAM';
            N_bps = 6;
            R_t = 5;
            R_b = 6;
        case 8
            modulation = '256-QAM';
            modulation0 = '256QAM';
            N_bps = 8;
            R_t = 3;
            R_b = 4;
        case 9
            modulation = '256-QAM';
            modulation0 = '256QAM';
            N_bps = 8;
            R_t = 5;
            R_b = 6;
        case 10
            modulation = '1024-QAM';
            modulation0 = '1024QAM';
            N_bps = 10;
            R_t = 3;
            R_b = 4;
        case 11
            modulation = '1024-QAM';
            modulation0 = '1024QAM';
            N_bps = 10;
            R_t = 5;
            R_b = 6;
        otherwise
            error('Unkown MCS index %f', MCS_index);
    end
    
    mcs.MCS_index = MCS_index;
    mcs.modulation = modulation;
    mcs.modulation0 = modulation0;
    mcs.N_bps = N_bps;
    mcs.R_t = R_t;
    mcs.R_b = R_b;
end
