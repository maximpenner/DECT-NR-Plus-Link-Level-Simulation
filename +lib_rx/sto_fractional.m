function [antenna_streams_mapped_rev_derot, sto_fractional] = sync_STO_fractional(antenna_streams_mapped_rev, physical_resource_mapping_STF_cell, N_RX, oversampling)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;

    % we need the size of the packet
    [N_b_DFT, N_PACKET_symb] = size(cell2mat(antenna_streams_mapped_rev(1)));

    % first get the position of all stf subcarriers and their values
    STF_l_matlab    = cell2mat(physical_resource_mapping_STF_cell(2)) + MATLAB_INDEX_SHIFT;
    STF_k_i_matlab  = lib_util.index_conversion_TS_matlab(N_b_DFT, cell2mat(physical_resource_mapping_STF_cell(1)));
    STF_values      = cell2mat(physical_resource_mapping_STF_cell(3));

    % our goal is to extract the angle between all sent and received STF cells
    ls_stf = zeros(numel(STF_values), N_RX);

    % for each rx antenna
    for i=1:1:N_RX

        % received f domain samples at antenna i
        transmit_streams_rev_i = cell2mat(antenna_streams_mapped_rev(i));

        % get angle between sent symbols and received symbols
        ls_stf(:, i) = STF_values .* conj(transmit_streams_rev_i(STF_k_i_matlab,STF_l_matlab));
    end

    %% extract the phase rotation between neighbouring subcarriers (this cancels out channel effects)
    ls_stf_s2s = ls_stf(1:end-1, :) .* conj(ls_stf(2:end, :));

    assert(mod(size(ls_stf_s2s, 1), 2) == 1);

    %% there are four subcarriers spacing between two occupied STF cells, expect to the center cells as can be seen in Figure 4.5-2.

    % what is the index where the spacing is not four subcarriers, but eight?
    center_idx = numel(STF_values) / 2;

    % it is the same complex value, but with only halve the angle
    ls_stf_s2s(center_idx, :) = abs(ls_stf_s2s(center_idx, :)) .* exp(1i * angle(ls_stf_s2s(center_idx, :)) / 2);

    %% average across all antennas (implicit weighting by magnitude)
    ls_stf_s2s_avg = sum(ls_stf_s2s, 'all');

    % subcarrier spacing between two occupied subcarrier is four
    phase_rotation = angle(ls_stf_s2s_avg) / 4;

    % convert phase rotation to absolute sample offset
    sto_fractional = phase_rotation / 2 / pi;
    sto_fractional = sto_fractional / (1 / (N_b_DFT * oversampling));

    %% create a vector to derotate signals of all antennas
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
        transmit_streams_rev_i_derot = transmit_streams_rev_i .* conj(derotation_mat);

        % write back
        antenna_streams_mapped_rev_derot(i) = {transmit_streams_rev_i_derot};
    end
end
