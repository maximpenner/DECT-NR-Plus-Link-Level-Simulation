function [y_PxC_ss] = Spatial_Multiplexing(x_PxC, N_SS)

    assert(N_SS > 1);
        
    % 6.3.2
    M_symb = numel(x_PxC);

    assert(mod(M_symb,N_SS) == 0);

    M_stream_symb = M_symb/N_SS;

    % de-interleave into spatial streams
    y_PxC_ss = cell(N_SS,1);
    for i=1:1:N_SS

        x_s = x_PxC(i:N_SS:end);

        assert(numel(x_s) == M_stream_symb);

        y_PxC_ss(i) = {x_s};
    end
end
