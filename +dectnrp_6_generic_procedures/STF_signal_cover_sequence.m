function [samples_antenna_tx_with_cover_sequence] = STF_signal_cover_sequence(samples_antenna_tx, ...
                                                                              u, ...
                                                                              b_times_oversampling)

    % lookup cover sequence
    if u==1
        c_u = [1; -1; 1; 1; -1; -1; -1];
        %c_u = ones(7,1);
    else
        c_u = [1; -1; 1; 1; -1; -1; -1; -1; -1];
        %c_u = ones(9,1);
    end

    % every pattern has a length of 16 samples, on top we put the oversampling
    c_u_16_os = repelem(c_u, 16*b_times_oversampling);

    % size of STF
    N_TX = size(samples_antenna_tx, 2);
    n_STF_samples_os = numel(c_u_16_os);

    % preallocate
    samples_antenna_tx_with_cover_sequence = samples_antenna_tx;

    for i=1:1:N_TX
        samples_antenna_tx_with_cover_sequence(1:n_STF_samples_os, i) = samples_antenna_tx(1:n_STF_samples_os, i).*c_u_16_os;
    end
end
