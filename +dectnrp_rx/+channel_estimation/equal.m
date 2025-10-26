function [weights_ts] = equal(N_b_DFT, ...
                              N_PACKET_symb, ...
                              N_p_used, ...
                              x_p, ...
                              y_p)
    % find indexes of the closest DRS pilots for every subcarrier
    idx_closest_drs = dectnrp_rx.channel_estimation.find_closest_drs(N_b_DFT, ...
                                                                      N_PACKET_symb, ...
                                                                      N_p_used, ...
                                                                      x_p, ...
                                                                      y_p);

    % return value
    weights_ts = zeros(N_b_DFT, N_PACKET_symb, numel(x_p));

    % go over each subcarrier ...
    for y_s = 1:1:N_b_DFT

        % ... in each OFDM symbol
        for x_s = 1:1:N_PACKET_symb
            
            % closest DRS pilots for this subcarrier
            idx_closest_drs_ = squeeze(idx_closest_drs(y_s, x_s, :));
            
            % write values at correct indexes
            weights_ts(y_s, x_s, idx_closest_drs_) = 1/N_p_used;
        end
    end

    assert(size(weights_ts, 1) == N_b_DFT);
    assert(size(weights_ts, 2) == N_PACKET_symb);
    assert(size(weights_ts, 3) == numel(x_p));
end
