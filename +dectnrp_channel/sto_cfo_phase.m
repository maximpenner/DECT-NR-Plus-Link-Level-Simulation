function [samples_antenna] = sto_cfo_phase(samples_antenna, sto_integer, sto_fractional, cfo, err_phase)

    N_TX = size(samples_antenna, 2);

    assert(sto_integer >= 0);

    % add integer sto
    samples_antenna = [zeros(sto_integer, N_TX); samples_antenna; zeros(sto_integer, N_TX)];

    % add fractional sto
    if sto_fractional ~= 0

        time_base = 0:size(samples_antenna, 1)-1;
        time_base = time_base';

        % delayseq() requires radar toolbox, interp1 is default functionality
        time_base_interpolation = time_base + sto_fractional;
        for i=1:1:N_TX
            samples_antenna(:, i) = interp1(time_base, samples_antenna(:, i), time_base_interpolation, 'spline', 'extrap');
        end
    end

    % add cfo
    if cfo ~= 0
        time_base = 0:size(samples_antenna, 1)-1;
        time_base = time_base';

        for i=1:1:N_TX
            samples_antenna(:, i) = samples_antenna(:, i) .* exp(1i*2*pi*cfo*time_base);
        end
    end

    % add error phase
    if err_phase ~= 0
        samples_antenna = samples_antenna*exp(1i*err_phase);
    end
end
