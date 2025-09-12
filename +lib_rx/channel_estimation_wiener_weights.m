function [wiener_weights_all_tx] = channel_estimation_wiener_weights(physical_resource_mapping_DRS_cell, ...
                                                                     N_closest_DRS_pilots, ...
                                                                     N_b_DFT, ...
                                                                     N_PACKET_symb, ...
                                                                     N_b_CP, ...
                                                                     samp_rate, ...
                                                                     noise_estim, ...
                                                                     f_d_hertz, ...
                                                                     tau_rms_sec)
    % frame parameters
    T_s_sec = 1/samp_rate;                      % sample time in seconds
    T_symb_sec = T_s_sec * (N_b_DFT + N_b_CP);  % length of symbol (CP included) in seconds
    f_subc_spacing_hertz = samp_rate/N_b_DFT;   % subcarrier spacing in Hertz
    
    % how many transmit antennas do we have?
    N_TX = size(physical_resource_mapping_DRS_cell,1);
    
    % how many pilots do we have for each tx antenna?
    N_p = numel(cell2mat(physical_resource_mapping_DRS_cell(1, end)));
    
    % how many of those pilots do we actually use?
    N_p_used = min(N_p, N_closest_DRS_pilots);
    
    % output, we save wiener filters for each transmit antenna
    wiener_weights_all_tx = cell(N_TX,1);
    
    %% go over each tx antenna as each tx antenna uses different pilots
    for tx_idx = 1:1:N_TX
        
        % First create an empty container for the wiener weights of this tx antenna.
        % If n_closest_pilots < N_p this matrix will contain zeros.
        wiener_weights_this_tx = zeros(N_b_DFT, N_PACKET_symb, N_p);
        
        % get linear indices of all pilots for this tx antenna
        linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(tx_idx, end));
        linear_indices_matlab = reshape(linear_indices_matlab,[],1);

        % convert linear indices to x and y coordinates
        [y_p, x_p] = ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab);

        % go over each subcarrier ...
        for i = 1:1:N_b_DFT

            % ... in each OFDM symbol
            for j = 1:1:N_PACKET_symb
                
                % get distances in tf lattice
                delta_t = x_p - j;
                delta_f = y_p - i;                
                
                % determine distances to all pilots
                distances = sqrt(delta_t.^2 + delta_f.^2);

                % find indices of closest pilots
                [~,idx_pilots_used] = mink(distances, N_closest_DRS_pilots);
                
                assert(numel(idx_pilots_used) == N_p_used);
                
                %% next we need to determine the correlation matrix between all used pilots
                R_pp = zeros(N_p_used,N_p_used);
                for ii = 1:1:N_p_used
                    
                    % extract position of this pilot
                    y_p_ref = y_p(idx_pilots_used(ii));
                    x_p_ref = x_p(idx_pilots_used(ii));

                    % what is the correlation to each used pilot
                    for jj = 1:1:N_p_used

                        % get relative difference
                        delta_t = x_p(idx_pilots_used(jj)) - x_p_ref;
                        delta_f = y_p(idx_pilots_used(jj)) - y_p_ref;

                        % set correlation values
                        corr_t = r_t(f_d_hertz, delta_t*T_symb_sec);
                        corr_f = r_f(tau_rms_sec, delta_f*f_subc_spacing_hertz);
                        corr_total = corr_t*corr_f;

                        % write to matrix
                        R_pp(jj,ii) = corr_total;
                    end
                end

                % we also need the noise on top of the estimation
                R_pp = R_pp + noise_estim*eye(N_p_used);
                
                %% next we need the correlation between this subcarrier and each pilot used
                R_dp = zeros(N_p_used, 1);
                for ii = 1:1:N_p_used

                    % extract position of this pilot
                    y_p_ = y_p(idx_pilots_used(ii));
                    x_p_ = x_p(idx_pilots_used(ii));

                    % get relative difference
                    delta_t = x_p_ - j;
                    delta_f = y_p_ - i;

                    % set correlation values
                    corr_t = r_t(f_d_hertz, delta_t*T_symb_sec);
                    corr_f = r_f(tau_rms_sec, delta_f*f_subc_spacing_hertz);
                    corr_total = corr_t*corr_f;

                    % write to matrix
                    R_dp(ii,1) =  corr_total;
                end
                
                %% solve wiener-hopf equation, weights are the wiener filter coefficients
                weights = R_pp\R_dp;

                % normalize, each channel coefficient has an expectation value of 1
                weights = weights/sum(weights);

                assert(abs(1 - sum(weights)) <= 1e-8);

                % write back
                wiener_weights_this_tx(i,j,idx_pilots_used) = weights;
            end
        end
        
        % save output of this tx antenna
        wiener_weights_all_tx(tx_idx) = {wiener_weights_this_tx};
    end 
end

% source: page 28 in f_subc_spacing https://publik.tuwien.ac.at/files/PubDat_204518.pdf
function correlation = r_t(nu_max_hertz, delta_t)
    correlation = besselj(0, 2 * pi * nu_max_hertz * delta_t);
end

% source: page 28 in f_subc_spacing https://publik.tuwien.ac.at/files/PubDat_204518.pdf
function correlation = r_f(tau_rms_sec, delta_f)
    correlation = 1/(1 + 1i * 2 * pi * tau_rms_sec * delta_f);
end
