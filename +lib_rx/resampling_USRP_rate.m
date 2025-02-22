function [USRP_samp_rate, L, M] = resampling_USRP_rate(DECT_samp_rate, oversampling, USRP_samp_rate_min_multiple)

    % constraints:
    %
    %   1) USRP_samp_rate > DECT_samp_rate_oversampled * USRP_samp_rate_min_multiple
    %
    %   2) L > M
    %

    DECT_samp_rate_oversampled = DECT_samp_rate*oversampling;

    % https://kb.ettus.com/USRP_N300/N310/N320/N321_Getting_Started_Guide
    USRP_master_clock   = 122.88e6;
    decimation_factor   = [2,4,6,8,10,12,14,16,18,20,30,32,64,100,128,200,256,512,1024];
    USRP_samp_rates     = USRP_master_clock./decimation_factor;

    % find sample rate above DECT_samp_rate_oversampled
    USRP_samp_rate = USRP_samp_rates(1);
    for i=2:1:numel(decimation_factor)

        if USRP_samp_rates(i) < DECT_samp_rate_oversampled * USRP_samp_rate_min_multiple
            break;
        else
            USRP_samp_rate = USRP_samp_rates(i);
        end
    end

    % sanity check
    if USRP_samp_rate < DECT_samp_rate_oversampled * USRP_samp_rate_min_multiple
        error("USRP samp rate smaller than DECT-2020 NR samp rate.");
    end

    % ratio between the sampling rates
    [L,M] = rat(USRP_samp_rate / DECT_samp_rate_oversampled);

    % sanity check
    % for b=1,2,4,8 and 16 is should be L=10 and M=9
    % for b=12 it should be L=40 and M=27
    if L == 10*USRP_samp_rate_min_multiple

        if M ~= 9
            error('Incorrect resampling coefficients.');
        end
    elseif L == 40*USRP_samp_rate_min_multiple

        if M ~= 27
            error('Incorrect resampling coefficients.');
        end
    else
            error('Incorrect resampling coefficients.');
    end
end