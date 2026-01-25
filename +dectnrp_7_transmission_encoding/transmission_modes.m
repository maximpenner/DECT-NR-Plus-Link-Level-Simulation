function [tm_mode] = transmission_modes(tm_mode_0_to_11)

    switch tm_mode_0_to_11
        case 0
            N_eff_TX = 1;
            N_SS = 1;
            CL = false;
            N_TS = 1;
            N_TX = 1;
            
        case 1
            N_eff_TX = 2;
            N_SS = 1;
            CL = false;
            N_TS = 2;
            N_TX = 2;
        case 2
            N_eff_TX = 2;
            N_SS = 2;
            CL = false;
            N_TS = 2;
            N_TX = 2;
        case 3
            N_eff_TX = 1;
            N_SS = 1;
            CL = true;
            N_TS = 1;
            N_TX = 2;
        case 4
            N_eff_TX = 2;
            N_SS = 2;
            CL = true;
            N_TS = 2;
            N_TX = 2;
            
        case 5
            N_eff_TX = 4;
            N_SS = 1;
            CL = false;
            N_TS = 4;
            N_TX = 4;
        case 6
            N_eff_TX = 4;
            N_SS = 4;
            CL = false;
            N_TS = 4;
            N_TX = 4;
        case 7
            N_eff_TX = 1;
            N_SS = 1;
            CL = true;
            N_TS = 1;
            N_TX = 4;
        case 8
            N_eff_TX = 2;
            N_SS = 2;
            CL = true;
            N_TS = 2;
            N_TX = 4;
        case 9
            N_eff_TX = 4;
            N_SS = 4;
            CL = true;
            N_TS = 4;
            N_TX = 4;
            
        case 10
            N_eff_TX = 8;
            N_SS = 1;
            CL = false;
            N_TS = 8;
            N_TX = 8;
        case 11
            N_eff_TX = 8;
            N_SS = 8;
            CL = false;
            N_TS = 8;
            N_TX = 8;
        otherwise
            error('Unknown transmission mode %f.', tm_mode_0_to_11);
    end
    
    tm_mode.mode_0_to_11 = tm_mode_0_to_11;
    tm_mode.N_eff_TX = N_eff_TX;
    tm_mode.N_SS = N_SS;
    tm_mode.CL = CL;
    tm_mode.N_TS = N_TS;
    tm_mode.N_TX = N_TX;    
end
