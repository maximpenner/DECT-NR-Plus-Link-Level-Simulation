function [samples_out] = err_phase_cfo_sto(samples_antenna, err_phase, cfo, sto)

    samples_out = samples_antenna;

    % add sto
    if sto >= 0
        samples_out = [zeros(sto,1); samples_out; zeros(sto,1)];
    else
        error('STO must be positive or zero, is %d.', sto);
    end 

    % add cfo
    if cfo ~= 0
        time_base = 0:numel(samples_out)-1;
        time_base = time_base';
        samples_out = samples_out.*exp(1i*2*pi*cfo*time_base);
    end

    % add error phase
    if err_phase ~= 0
        samples_out = exp(1i*err_phase)*samples_out;
    end
end