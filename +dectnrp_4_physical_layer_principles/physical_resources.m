function [k_b_OCC] = physical_resources(N_b_OCC)

    % occupied subcarriers
    k_b_OCC = -N_b_OCC/2 : 1 : N_b_OCC/2;

    % remove DC
    k_b_OCC(N_b_OCC/2+1) = [];
end
