function [weights_ts] = wiener(N_b_DFT, ...
                               N_PACKET_symb, ...
                               N_b_CP, ...
                               N_p_used, ...
                               x_p, ...
                               y_p, ...
                               samp_rate, ...
                               noise_estim, ...
                               f_d_hertz, ...
                               tau_rms_sec)
    % packet parameters
    T_s_sec = 1/samp_rate;
    T_symb_sec = T_s_sec * (N_b_DFT + N_b_CP);
    f_subc_spacing_hertz = samp_rate/N_b_DFT;

    % return value
    weights_ts = zeros(N_b_DFT, N_PACKET_symb, numel(x_p));

    % go over each subcarrier ...
    for y_s = 1:1:N_b_DFT

        % ... in each OFDM symbol
        for x_s = 1:1:N_PACKET_symb
            
            % get distances in tf lattice
            delta_t = x_p - x_s;
            delta_f = y_p - y_s;                
            
            % determine distances to all pilots
            distances = sqrt(delta_t.^2 + delta_f.^2);

            % find indices of closest pilots
            [~,idx_pilots_used] = mink(distances, N_p_used);
            
            assert(numel(idx_pilots_used) == N_p_used);
            
            %% next we need to determine the correlation matrix between all pilots used
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
                delta_t = x_p_ - x_s;
                delta_f = y_p_ - y_s;

                % set correlation values
                corr_t = r_t(f_d_hertz, delta_t*T_symb_sec);
                corr_f = r_f(tau_rms_sec, delta_f*f_subc_spacing_hertz);
                corr_total = corr_t*corr_f;

                % write to matrix
                R_dp(ii,1) =  corr_total;
            end
            
            %% solve Wiener-Hopf equation, weights are the wiener filter coefficients
            weights = R_pp\R_dp;

            % normalize, each channel coefficient has an expectation value of 1
            weights = weights/sum(weights);

            assert(abs(1 - sum(weights)) <= 1e-8);

            weights_ts(y_s, x_s, idx_pilots_used) = weights;
        end
    end

    assert(size(weights_ts, 1) == N_b_DFT);
    assert(size(weights_ts, 2) == N_PACKET_symb);
    assert(size(weights_ts, 3) == numel(x_p));
end

function correlation = r_t(nu_max_hertz, delta_t)
    % source: page 28 in f_subc_spacing https://publik.tuwien.ac.at/files/PubDat_204518.pdf
    correlation = besselj(0, 2 * pi * nu_max_hertz * delta_t);
end

function correlation = r_f(tau_rms_sec, delta_f)
    % source: page 28 in f_subc_spacing https://publik.tuwien.ac.at/files/PubDat_204518.pdf
    correlation = 1/(1 + 1i * 2 * pi * tau_rms_sec * delta_f);
end
