function [samples_antenna_stf_at_coarse_peak_cfo_fractional, cfo_fractional_report] = cfo_fractional(samples_antenna_stf_at_coarse_peak, ...
                                                                                                     n_STF_template, ...
                                                                                                     u)
    % SOURCES
    % Schmidl Cox: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=650240
    % BLUE: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=752907
    % wiki: https://en.wikipedia.org/wiki/Carrier_frequency_offset
    
    % how often does the pattern repeat?
    switch u
        case 1
            L = 7;
        case {2,4,8}
            L = 9;
    end

    assert(mod(n_STF_template, L*16) == 0, 'oversampling not an integer value');

    % n_STF_template is the oversampled length of the STF, by dividing with L*16 we get the oversampling itself
    oversampling = n_STF_template / (L*16);

    % revert STF cover sequence by applying it
    samples_antenna_stf_at_coarse_peak = dectnrp_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_stf_at_coarse_peak, u, oversampling);
    
    % number of samples in a single pattern repetition
    n_STF_pattern = n_STF_template/L;
    
    % get the size of the input packet
    [n_samples_antenna_sto, N_RX] = size(samples_antenna_stf_at_coarse_peak);

    % output container
    samples_antenna_stf_at_coarse_peak_cfo_fractional = zeros(size(samples_antenna_stf_at_coarse_peak));
    
    % time base for CFO correction
    time_base = 0:1:(n_samples_antenna_sto-1);
    time_base = time_base';

    % equation (5) from Schmidl Cox: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=650240
    P_total = 0;
    
    % process each antenna
    for i=1:1:N_RX
        
        % extract the preamble as we have received it, perfect symbol synchronization is assumed
        STF_found = samples_antenna_stf_at_coarse_peak(1:n_STF_template, i);

        % SCHMIDL COX
        % Schmidl Cox assumes one OFDM symbol that consists of 2 equal halves.
        % DECT-2020 has 1.75 (u=1) or 2.25 (u=2,4,8) OFDM symbols with 7 or 9 repetitions.
        % Our capture rate is therefore doubled in comparison to Schmidl Cox, we can capture a CFO of +/- two subcarrier spacings.
        % Schmidl Cox can capture only +/- one subcarrier spacing.

        % reorder to repeated patterns
        STF_found_mat = reshape(STF_found, n_STF_pattern, L);
        
        % equation (5)
        P = STF_found_mat(:,2:end).*conj(STF_found_mat(:,1:end-1));

        % summation equation (5)
        P = sum(sum(P));
        
        % estimate the cfo across all antennas
        P_total = P_total + P;
    end

    cfo_fractional_report = cfo_from_p_SchmidlCox(P_total, n_STF_pattern);
    
    % derotate carrier frequency offset for each antenna with one common CFO estimation
    for i=1:1:N_RX
        samples_antenna_stf_at_coarse_peak_cfo_fractional(:,i) = samples_antenna_stf_at_coarse_peak(:,i).*exp(1i*2*pi*(-cfo_fractional_report)*time_base);
    end
end

function cfo_report = cfo_from_p_SchmidlCox(P, n_STF_pattern)
    
    % by how much does the fractional STO rotate over one pattern normalized to 2*pi?
    cfo_report = angle(P)/(2*pi);

    % this line can be used to see the rotation over one OFDM symbol without CP
    %cfo_report_one_symbol = cfo_report*4;

    % The CFO we need to use for correction at a sampling rate of 1.
    % Note that N_b_DFT_os = 4*n_STF_pattern.
    cfo_report = cfo_report/n_STF_pattern;
end
