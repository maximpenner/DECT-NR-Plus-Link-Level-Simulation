function [fir_coef_tx, fir_coef_rx] = resampling_filter_design( L, ...
                                                                M, ...
                                                                f_pass_norm, ...
                                                                f_stop_norm, ...
                                                                passband_ripple_dB, ...
                                                                stopband_attenuation_dB)
    % sanity checks
    if L <= M
        error('For resampling L must be larger than M.');
    end
    if f_pass_norm >= 0.5
        error('For critically sampled signals passband must end below 0.5 Hz.');
    end
    if f_pass_norm >= f_stop_norm
        error('f_pass_norm must be smaller than f_stop_norm.');
    end

    % Image filter and ani aliasing filter merge into one.
    % Since L > M, we use L at both TX and RX for determining the passband, as this creates a narrower filter.

    % TX resampling filter design

    f_pass = f_pass_norm * 1/L;
    f_stop = f_stop_norm * 1/L;

    % generate a filter with Matlab function
    fir_tx = designfilt('lowpassfir', ...
                        'PassbandFrequency',f_pass, ...
                        'StopbandFrequency',f_stop, ...
                        'PassbandRipple',passband_ripple_dB, ...
                        'StopbandAttenuation',stopband_attenuation_dB, ...
                        'DesignMethod','kaiserwin', ...
                        'SampleRate',1);

    % extract filter coefficients
    fir_coef_tx_reference = L*tf(fir_tx);

    % generate filter coefficients with custom function
    fir_coef_tx = lib_rx.resampling_kaiser( fir_tx.PassbandFrequency,...
                                            fir_tx.StopbandFrequency,...
                                            fir_tx.PassbandRipple,...
                                            fir_tx.StopbandAttenuation,...
                                            fir_tx.SampleRate,...
                                            true);

    fir_coef_tx = L*fir_coef_tx;

    % RX resampling filter design

    f_pass = f_pass_norm * 1/L;
    f_stop = f_stop_norm * 1/L;

    % generate a filter with Matlab function
    fir_rx = designfilt('lowpassfir', ...
                        'PassbandFrequency',f_pass, ...
                        'StopbandFrequency',f_stop, ...
                        'PassbandRipple',passband_ripple_dB, ...
                        'StopbandAttenuation',stopband_attenuation_dB, ...
                        'DesignMethod','kaiserwin', ...
                        'SampleRate',1);

    % extract filter coefficients
    fir_coef_rx_reference = M*tf(fir_rx);

    % generate filter coefficients with custom function
    fir_coef_rx = lib_rx.resampling_kaiser( fir_rx.PassbandFrequency,...
                                            fir_rx.StopbandFrequency,...
                                            fir_rx.PassbandRipple,...
                                            fir_rx.StopbandAttenuation,...
                                            fir_rx.SampleRate,...
                                            true);

    fir_coef_rx = M*fir_coef_rx;

    % debug
    if 1==0

        % TX filter

        % plot filter with matlab function
        fvtool(fir_tx)
    
        % plot filters im time domain for comparison
        figure()
        clf();

        % reference filter matlab
        N = numel(fir_coef_tx_reference);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_tx_reference, 'LineWidth', 2.2)
        hold on

        % custom filter
        N = numel(fir_coef_tx);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_tx, 'r')
        legend('Matlab Reference', 'Custom Function')
        title('TX Resampling Filter')

        % RX filter

        % plot filter with matlab function
        fvtool(fir_rx)

        % plot filters im time domain for comparison
        figure()
        clf();

        % reference filter matlab
        N = numel(fir_coef_rx_reference);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_rx_reference, 'LineWidth', 2.2)
        hold on

        % custom filter
        N = numel(fir_coef_rx);
        n = 0:1:(N-1);
        plot(n-(N-1)/2, fir_coef_rx, 'r')
        legend('Matlab Reference', 'Custom Function')
        title('RX Resampling Filter')
    end
end

