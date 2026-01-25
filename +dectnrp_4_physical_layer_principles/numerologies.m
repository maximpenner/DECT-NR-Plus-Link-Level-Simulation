function [numerology] = numerologies(u, b)

    assert(ismember(u,[1, 2, 4, 8]));
    
    assert(ismember(b,[1, 2, 4, 8, 12, 16]));
    
    % Table 4.3-1
    delta_u_f = u*27000;
    T_u_symb = 1/delta_u_f * (64+8)/64; % (64+8)/64 = (128+16)/128 = (256+32)/256 = ...
    N_SLOT_u_symb = u*10;
    N_SLOT_u_subslot = u*2;
    
    switch u
        case 1
            GI_u = 4/9*T_u_symb;        % Figure 5.1-1
        case {2,4}
            GI_u = T_u_symb;            % Figure 5.1-2
        case 8
            GI_u = 2*T_u_symb;          % Figure 5.1-3
    end
    
    B_u_b_DFT = b*64*delta_u_f;
    T_u_b_s = 1/B_u_b_DFT;
    N_b_DFT = b*64;
    N_b_CP = b*8;
    N_b_OCC = b*56;
    B_u_b_TX = N_b_OCC*delta_u_f;
    
    % add everything to output
    numerology.delta_u_f        = delta_u_f;
    numerology.T_u_symb         = T_u_symb;
    numerology.N_SLOT_u_symb    = N_SLOT_u_symb;
    numerology.N_SLOT_u_subslot = N_SLOT_u_subslot;
    numerology.GI_u             = GI_u;
    numerology.B_u_b_DFT        = B_u_b_DFT;
    numerology.T_u_b_s          = T_u_b_s;
    numerology.N_b_DFT          = N_b_DFT;
    numerology.N_b_CP           = N_b_CP;
    numerology.N_b_OCC          = N_b_OCC;
    numerology.B_u_b_TX         = B_u_b_TX;

    % we add a few more useful parameters

    % also calculate the guards
    numerology.n_guards_top = (N_b_DFT - N_b_OCC)/2 - 1;
    numerology.n_guards_bottom = numerology.n_guards_top + 1;

    assert(numerology.n_guards_bottom + N_b_OCC + numerology.n_guards_top + 1 == N_b_DFT);
end
