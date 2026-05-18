function [antenna_streams_mapped_rev_derot, sto_residual] = sto_residual(antenna_streams_mapped_rev, ...
                                                                         physical_resource_mapping_DRS_cell, ...
                                                                         physical_resource_mapping_STF_cell, ...
                                                                         N_RX, ...
                                                                         N_eff_TX, ...
                                                                         oversampling)
    % After FFT a remaining STO appears as a phase slope across the subcarriers,
    % this rotation can be estimated by neighboring pilots where H_k1 * H_k2 ≈ |H_k1|^2.
    % The equation is from https://ieeexplore.ieee.org/document/4017715.

    % we need the dimensions of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    %% STF preparation (optional/ already done in sto_fractional)

    % The STF could also be used for this algorithm, since it is scheduled
    % in a similar manor (every forth subcarrier) as the pilots, however,
    % sto correction based on the STF after FFT is already done in sto_fractional.

    % For each RX antenna, the same STF is received.
    % Get the STF indices and the sent STF complex values.
    STF_linear_indices_matlab = dectnrp_util.index_conversion_TS_matlab(N_b_DFT, physical_resource_mapping_STF_cell{1});
    STF_values                = physical_resource_mapping_STF_cell{3};

    %% DRS preparation

    DRS_subcarrier_indices_i = cell2mat(physical_resource_mapping_DRS_cell(1,1));

    % Keep the signed spacing between neighboring occupied DRS subcarriers.
    % This sign is required for a correct STO sign in the phase-slope estimate.
    dk = DRS_subcarrier_indices_i(2,1) - DRS_subcarrier_indices_i(1,1);

    %% output preparation

    % init phase comparison
    a = 0;

    %% the residual STO is corrected for each RX antenna individually
    for i=1:1:N_RX

        % extract received frequency domain samples for RX antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % STF values already in use at STO fractional

        %% extract received complex STF values (sent values were extracted above)

        % received
        y_STF = transmit_streams_rev_i(STF_linear_indices_matlab);

        s = 1:2:numel(y_STF);

        a = a + sum((y_STF(s).* conj(STF_values(s))).*conj(y_STF(s+1).* conj(STF_values(s+1)))); 

        for j=1:1:N_eff_TX

            %% extract sent and received DRS values

            % linear indices of drs values
            DRS_linear_indices_matlab = cell2mat(physical_resource_mapping_DRS_cell(j,4));

            % transmitted drs
            DRS_values = cell2mat(physical_resource_mapping_DRS_cell(j,3));

            % received drs
            y_DRS = transmit_streams_rev_i(DRS_linear_indices_matlab);

            % In Eq. (19), accumulate (k1) * conj(k2) with dk = k2 - k1.
            for l=1:size(y_DRS,2)
                for s = 1:2:numel(DRS_values)
                    a = a + (y_DRS(s,l) * conj(DRS_values(s))) * conj(y_DRS(s+1,l) * conj(DRS_values(s+1)));
                end
            end
        end
    end

    % Eq. (19): theta_hat = N/(2*pi*dk) * angle(sum(...))
    % Thus phase slope across subcarriers is: phase_rotation = 2*pi*theta_hat/N = angle(sum(...))/dk
    phase_rotation = angle(a) / dk;

    % Report in time-domain sample units (same convention as sto_fractional).
    sto_residual = phase_rotation / 2 / pi;
    sto_residual = sto_residual / (1 / (N_b_DFT * oversampling));

    %% create a vector to derotate signals of all antennas
    % Row order in mapped streams follows TS subcarrier indices from +N/2-1 down to -N/2.
    % Keep this order to apply the phase-slope correction to the correct bins.
    % k_vec_TS_order = (N_b_DFT/2-1) : -1 : (-N_b_DFT/2);
    % derotation_vec = exp(1i * k_vec_TS_order.' * phase_rotation);
    % derotation_mat = repmat(derotation_vec, 1, N_PACKET_symb);
    derotation_vec = -N_b_DFT/2 : 1 : (N_b_DFT/2-1);
    derotation_vec = derotation_vec';
    derotation_vec = exp(1i*derotation_vec*(-phase_rotation));
    derotation_mat = repmat(derotation_vec, 1, N_PACKET_symb);

    %% create a copy of the input signal and then derotate that copy
    antenna_streams_mapped_rev_derot = antenna_streams_mapped_rev;

    % derotate signal of each antenna
    for i=1:1:N_RX
        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % derotate
        transmit_streams_rev_i_derot = transmit_streams_rev_i .* derotation_mat;

        % write back
        antenna_streams_mapped_rev_derot(i) = {transmit_streams_rev_i_derot};
    end
end
