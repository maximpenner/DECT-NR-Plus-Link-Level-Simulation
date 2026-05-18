function [codebook_index_feedback,N_TS_feedback] = MIMO_feedback(antenna_streams_mapped_rev, ...
                                                                 physical_resource_mapping_DRS_cell, ...
                                                                 ch_estim, ...
                                                                 N_RX, ...
                                                                 N_eff_TX, ...
                                                                 noise_estim)
    % Calculation of RI, PMI and CQI is mostly based on Singular Value
    % Decomposition
    % Good description: https://www.sharetechnote.com/html/BasicProcedure_LTE_MIMO.html
    % Basically the Channel Matrix can be divided into three matrices
    %   V = precoding matrix
    %   U = receiver matrix
    %   S = diagonal matrix, which can be used to calculate the Rank Indicator = Number of Nonzero Values in S
    
    % [U, S, V] = pagesvd(H_est, "econ");
    
    % Why is this not used? Calculation of the SVD has to be done for every
    % subcarrier, this is very costly
    
    % Practically the RI and the PMI are selected according to a maximum calculated capacity
    % E.G. https://ieeexplore.ieee.org/document/5670402
    
    % In https://ieeexplore.ieee.org/document/6691823 and https://ieeexplore.ieee.org/document/8059870
    % a good calculation is shown based on LTE by using the channel correlation matrix which can be
    % averaged over the whole bandwidth, for this we can use the pilots (so we use the code from
    % +/dectnrp_rx/channel_estimation/estimate.m)
    %
    % Also we need a good noise estimation, otherwise the capacity calculation won't work.

    %% calculate correlation matrix (or use ch_estim instead)

    assert(N_RX >= N_eff_TX, "number of transmit and receive antennas must be equal for MIMO!");

    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    [ls_k ,ls_m] = size(cell2mat(physical_resource_mapping_DRS_cell(1,1)));

    % channel containers
    H_est = zeros(N_RX, N_eff_TX, ls_k, ls_m);
    ch_estim_size = size(ch_estim{1});
    ch_estim_size = [ch_estim_size , numel(ch_estim)];
    ch_estim_full = zeros(ch_estim_size);

    % for each RX antenna
    for i=1:1:N_RX

        ch_estim_i = ch_estim{i};

        ch_estim_full(:,:,:,i) = ch_estim_i;

        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % for each transmit stream
        for j=1:1:N_eff_TX
            %% extract send values
            % linear indices of DRS
            linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted DRS
            values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received DRS
            y = transmit_streams_rev_i(ind2sub([N_b_DFT N_PACKET_symb], linear_indices_matlab));

            %% least squares estimation at the pilot positions
            ls = y./repmat(values,1, size(y,2));

            H_est(i,j,:,:) = ls;
        end
    end

    R_H = zeros(N_RX, N_eff_TX, ls_m);
    
    % calculate channel correlation matrix
    for m = 1:ls_m
        for k=1:ls_k
            R_H(:,:,m) = R_H(:,:,m) + H_est(:,:,k,m)' * H_est(:,:,k,m);
        end
        R_H(:,:,m) = R_H(:,:,m)./ ls_k;
    end

    % permute to N_RX x N_TS x N_b_DFT x N_PACKET_symb
    ch_estim_full = permute(ch_estim_full, [4 3 1 2]);
    ch_estim_full = reshape(ch_estim_full, N_RX, N_eff_TX, []);
    
    %% look for the possible candidates for spatial streams (Rank Indicator) and codebook index
    switch N_eff_TX
        case 1
            error("MIMO not possible with one N_TX");
        case {2, 4}
            N_TS_candidates = 2.^(0:log2(N_eff_TX));
        case 8     
            N_TS_feedback = 8;
            codebook_index_feedback = 0;
    end
    
    % for 8x8 there is no closed loop mimo option
    if N_eff_TX < 8
        % maximum of transmit stream is dependent on number of receive antennas
        while (N_TS_candidates(end) > N_RX)
            N_TS_candidates(end) = [];
        end
    
        % look for all possible codebook indexes we could use
        codebook_index_candidates = cell(1, numel(N_TS_candidates));
        for i = 1:numel(N_TS_candidates)
            switch N_TS_candidates(i)
                case 1
                    switch N_eff_TX
                        case 1
                            error("MIMO not possible")
                        case 2
                            codebook_index_candidates{i} = 0:1:5;
                        case 4
                            codebook_index_candidates{i} = 0:1:27;
                    end
    
                case 2
                    switch N_eff_TX
                        case 1
                             error("MIMO not possible")
                        case 2
                            codebook_index_candidates{i} = 0:1:2;
                        case 4
                            codebook_index_candidates{i} = 0:1:21;
                    end
                    
                case 4
                    if N_eff_TX ~= 4
                        error(" four transmit streams requires 4 antennas")
                    else
                        codebook_index_candidates{i} = 0:1:4;
                    end
            
                otherwise
                    error("number of transmit streams not possible")
            end
        end

        %% calculate the metric and pick the best candidate
        % calculate the MMSE capacity according to https://ieeexplore.ieee.org/document/5670402
        % calculate the capacity for every subband/ symbol channel
        % estimation
        % even the Guard Band subcarriers can be used (at least it seems to
        % work)
        % average the capacity over the whole bandwidth and find the RI/PMI
        % pair where the minimum is

        W_k_m_sum = zeros(1, numel(N_TS_candidates));
        codebook_index_candidates_N_TS = zeros(1, numel(N_TS_candidates));

        for i=1:numel(N_TS_candidates)
            current_r = N_TS_candidates(i);
            codebook_index_candidates_i = codebook_index_candidates{i};
            c_full = zeros(1, numel(codebook_index_candidates_i));

            for j = 1:numel(codebook_index_candidates_i)
                current_W = dectnrp_6_generic_procedures.Beamforming_W(N_TS_candidates(i), N_eff_TX, codebook_index_candidates_i(j));

                % either the correlation matrix or ch_estim can be
                % used. ch_estim seems to work better (at 0 dB SNR a good
                % precoding matrix is selected)

                % for m = 1:numel(R_H(1,1,:))
                % 
                %     R_H_m = R_H(:,:,m);
                % 
                %     M = current_W' * R_H_m * current_W * 1/(noise_estim.^2) + eye(current_r);
                %     invM = inv(M);
                % 
                %     gamma = prod(diag(invM));
                % 
                %     I_k_m = sum(log2(real(gamma)));
                % 
                %     W_k_m_sum_j(j) = W_k_m_sum_j(j) + (I_k_m / numel(R_H(1,1,:)));
                % end

                for m = 1:numel(ch_estim_full(1,1,:))

                    ch_estim_full_m = ch_estim_full(:,:,m);

                    M = current_W' * (ch_estim_full_m' * ch_estim_full_m) * current_W * 1/(noise_estim.^2) + eye(current_r);
                    invM = inv(M);

                    gamma = prod(diag(invM));

                    c_m = sum(log2(real(gamma)));

                    c_full(j) = c_full(j) + (c_m / numel(ch_estim_full(1,1,:)));
                end
            end

            % find best candidate codebook index for the current candidate
            % N_TS (i)

            [~,min_index_W_k_m] = min(c_full);

            % fail safe (if multiple minimums, pick the first)
            if (numel(min_index_W_k_m) > 1)
                min_index_W_k_m_x = find(min_index_W_k_m);
                min_index_W_k_m_x = min_index_W_k_m_x(1);
                min_index_W_k_m(min_index_W_k_m_x+1:end) = [];
            end

            codebook_index_candidates_N_TS(i) = codebook_index_candidates_i(min_index_W_k_m);
            
            W_k_m_sum(i) = c_full(min_index_W_k_m);
        end

        % compare the result of candidate transmit streams

        [~,min_index_W] = min(W_k_m_sum);

        % Precoding Matrix Indicator
        codebook_index_feedback = codebook_index_candidates_N_TS(min_index_W);

        % Rank Indicator
        N_TS_feedback = N_TS_candidates(min_index_W);
    end
end