function [idx_closest_drs] = find_closest_drs(N_b_DFT, ...
                                              N_PACKET_symb, ...
                                              N_p_used, ...
                                              x_p, ...
                                              y_p)
    % return value
    idx_closest_drs = zeros(N_b_DFT, N_PACKET_symb, N_p_used);

    % go over each subcarrier ...
    for y_s = 1:1:N_b_DFT

        % ... in each OFDM symbol
        for x_s = 1:1:N_PACKET_symb
            
            % get distances to all pilots in tf lattice
            delta_t = x_p - x_s;
            delta_f = y_p - y_s;
            
            % determine distances to all pilots
            distances = sqrt(delta_t.^2 + delta_f.^2);

            % find indices of closest pilots
            [~, tmp] = mink(distances, N_p_used);

            idx_closest_drs(y_s, x_s, :) = tmp;
        end
    end

    assert(size(idx_closest_drs, 1) == N_b_DFT);
    assert(size(idx_closest_drs, 2) == N_PACKET_symb);
    assert(size(idx_closest_drs, 3) == N_p_used);
end
