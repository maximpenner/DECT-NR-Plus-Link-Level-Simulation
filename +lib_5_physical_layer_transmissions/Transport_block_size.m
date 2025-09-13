function [N_TB_bits] = Transport_block_size(tm_mode, mcs, N_PDC_subc, Z)

    % copy all relevant variables
    N_SS = tm_mode.N_SS;
    N_bps = mcs.N_bps;
    R_t = mcs.R_t;
    R_b = mcs.R_b;
    
    %% 5.3
    N_PDC_bits = floor(N_SS * N_PDC_subc * N_bps * R_t / R_b);
    
    L = 24;
    
    if N_PDC_bits <= 512
        M = 8;
    elseif N_PDC_bits <= 1024
        M = 16;
    elseif N_PDC_bits <= 2048
        M = 32;
    else
        M=64;
    end
    
    N_M = floor(N_PDC_bits/M) * M;
    
    if N_M <= Z
        N_TB_bits = N_M - L;
    else
        C = ceil((N_M-L)/Z);
        N_TB_bits = N_M - (C+1)*L;
    end
end
