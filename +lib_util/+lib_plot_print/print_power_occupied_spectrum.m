function [] = print_power_occupied_spectrum(antenna_streams_mapped, n_spectrum_occupied, samples_antenna_tx, n_STF_samples)

    % what is the power in time domain across all antennas?
    n_samples = size(samples_antenna_tx,1);
    power_across_all_antennas_t_domain = sum(sum(abs(samples_antenna_tx).^2))/n_samples;
    
    % what is the power in time domain across all antennas?
    power_stf_across_all_antennas_t_domain = sum(sum(abs(samples_antenna_tx(1:n_STF_samples,:)).^2))/n_STF_samples;    

    % measure the number of non-zero elements in all antenna streams
    antenna_streams_total_subcarriers = 0;
    occupied_subcarriers = 0;
    for i=1:1:numel(antenna_streams_mapped)
        temp = cell2mat(antenna_streams_mapped(i));
        antenna_streams_total_subcarriers = antenna_streams_total_subcarriers + numel(temp);
        occupied_subcarriers = occupied_subcarriers + sum(sum(temp ~= 0));
    end
    occupied_spectrum = occupied_subcarriers/antenna_streams_total_subcarriers;
    
    % show info
    disp('##### TX Packet ######');
    fprintf('Total power across all antennas in t domain:        %f\n', power_across_all_antennas_t_domain);
    fprintf('Total power of STF across all antennas in t domain: %f\n', power_stf_across_all_antennas_t_domain);
    fprintf('Occupied spectrum:                                  %f\n', occupied_spectrum);
    fprintf('Occupied spectrum scaled by n_spectrum_occupied:    %f\n\n', occupied_spectrum/n_spectrum_occupied);
end
