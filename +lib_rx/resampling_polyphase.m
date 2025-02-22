function samples_antenna_tx_or_rx_resampled = resampling_polyphase(samples_antenna_tx_or_rx, L, M, fir_coef_tx_or_rx)

    %% OPTION A: Matlab function
    samples_antenna_tx_or_rx_resampled_reference = resample(samples_antenna_tx_or_rx, L, M, fir_coef_tx_or_rx);

    %% Option B: custom function which can be implemented in C++

    % FIR filter setup

    n_samples_fir = numel(fir_coef_tx_or_rx);

    % sanity check
    if mod(n_samples_fir,2) == 0
        error('Number of input filter coefficients should be odd to have an integer delay.');
    end

    % what is the delay of the fir filter?
    n_delay = (n_samples_fir-1)/2;

    % make sure fir_coef_tx_or_rx has a length which is a multiple of L, append zeros if required
    n_fir_append_zeros = ceil(n_samples_fir/L)*L - n_samples_fir;
    fir_coef_tx_or_rx_ZEROS = [fir_coef_tx_or_rx, zeros(1,n_fir_append_zeros)];

    % now create subfilters
    n_subfilters_coeff = (n_samples_fir + n_fir_append_zeros)/L;
    subfilters = zeros(L, n_subfilters_coeff);
    for i=1:1:L
        subfilters(i,:) = fir_coef_tx_or_rx_ZEROS(i:L:end);
    end

    % Input processing

    [n_samples_input, N_RX] = size(samples_antenna_tx_or_rx);

    % what is the final number of output samples we require?
    % source: https://de.mathworks.com/help/signal/ref/resample.html#bumhz33-y
    n_samples_output = ceil(n_samples_input*L/M);

    % create output container
    samples_antenna_tx_or_rx_resampled = zeros(n_samples_output, N_RX);

    % we have to add some zeros at the beginning of the input samples as history
    n_samples_input_prepend_zeros = n_subfilters_coeff-1;

    % we also have to add some zeros at the end to make sure we have the correct number of output samples
    n_samples_input_append_zeros = ceil((n_delay+1 + (n_samples_output-1)*M)/L) - n_samples_input;

    % process each antenna signal
    for i=1:1:N_RX
    
        % add zeros to input samples
        samples_antenna_tx_or_rx_ZEROS = [  zeros(n_samples_input_prepend_zeros, 1); ...
                                            samples_antenna_tx_or_rx(:, i); ...
                                            zeros(n_samples_input_append_zeros,1)];
    
        % how apply each subfilter individually
        samples_upsampled = zeros(L, n_samples_input + n_samples_input_append_zeros);
        for j=1:1:L
    
            fir = subfilters(j,:);
            fir = fir';
    
            samples_upsampled(j,:) = conv(samples_antenna_tx_or_rx_ZEROS, fir, 'valid');
        end
    
        % now extract the downsampled samples, skip the filter delay
        samples_upsampled = reshape(samples_upsampled, [], 1);
        samples_downsampled = samples_upsampled(n_delay+1:M:end);
    
        % set final output
        samples_antenna_tx_or_rx_resampled(:, i) = samples_downsampled;
    
        %% sanity check
        err = max(abs(samples_antenna_tx_or_rx_resampled_reference(:,i) - samples_antenna_tx_or_rx_resampled(:, i)));
        if err > 1e-5
            error('Incorrect resampling results.');
        end
        
        %% DEBUG plot
        if 1==0
    
            % create time bases
            t_original  = 0:1:(numel(samples_antenna_tx_or_rx)-1);
            t_resampled = 0:1:(numel(samples_antenna_tx_or_rx_resampled_reference)-1);
            t_resampled = t_resampled/L*M;
        
            figure()
            clf()
        
            % plot real part of both signals
            subplot(4,1,1)
            plot(t_original, real(samples_antenna_tx_or_rx))
            hold on
            plot(t_resampled, real(samples_antenna_tx_or_rx_resampled_reference),'r')
            ylabel('REAL reference');
            title('Interpolation (blue original, red interpolated)')
            ylim([-3 3]);
        
            % plot imag part of both signals
            subplot(4,1,2)
            plot(t_original, imag(samples_antenna_tx_or_rx))
            hold on
            plot(t_resampled, imag(samples_antenna_tx_or_rx_resampled_reference),'r')
            ylabel('IMAG reference');
            title('Interpolation (blue original, red interpolated)')
            ylim([-3 3]);
    
            % plot real part of both signals
            subplot(4,1,3)
            plot(t_original, real(samples_antenna_tx_or_rx))
            hold on
            plot(t_resampled, real(samples_antenna_tx_or_rx_resampled),'r')
            ylabel('REAL');
            title('Interpolation (blue original, red interpolated)')
            ylim([-3 3]);
        
            % plot imag part of both signals
            subplot(4,1,4)
            plot(t_original, imag(samples_antenna_tx_or_rx))
            hold on
            plot(t_resampled, imag(samples_antenna_tx_or_rx_resampled),'r')
            ylabel('IMAG');
            title('Interpolation (blue original, red interpolated)')
            ylim([-3 3]);
        end
    end
end