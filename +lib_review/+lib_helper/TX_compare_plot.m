function [] = TX_compare_plot(  filename, ...
                                samples_antenna_tx_matlab, ...
                                samples_antenna_tx_matlab_resampled, ...
                                samples_antenna_tx_cpp, ...
                                dect_samp_rate, ...
                                hw_samp_rate)

    % for plotting, create two time axis (can be the same is no resampling is used)
    t_dect = (0:1:(numel(samples_antenna_tx_matlab(:,1)) - 1) ) * 1/dect_samp_rate;
    t_hw   = (0:1:(numel(samples_antenna_tx_matlab_resampled(:,1)) - 1) ) * 1/hw_samp_rate;

    % plot comparison
    screen_size_px = get(0,'screensize');
    screen_size_x = screen_size_px(3);
    screen_size_y = screen_size_px(4);
    x_margin = 50;
    y_margin = 50;
    for i=1:1:size(samples_antenna_tx_matlab, 2)

        fig = figure(i);

        % figure position
        if i <= 4
            fig_pos_x = x_margin + (screen_size_x - 2*x_margin) / 4 * (i-1);
            fig_pos_y = -y_margin;
            fig_pos = [fig_pos_x fig_pos_y];
        else
            fig_pos_x = x_margin + (screen_size_x - 2*x_margin) / 4 * (i-5);
            fig_pos_y = -(y_margin + screen_size_y / 2);
            fig_pos = [fig_pos_x fig_pos_y];
        end
        movegui(fig,fig_pos);
        clf()

        % compare real part
        subplot(5,1,1)
        plot(real(samples_antenna_tx_matlab_resampled(:,i)));
        hold on
        plot(real(samples_antenna_tx_cpp(:,i)), 'r');
        legend('Matlab', 'C++');
        grid on

        t = title(strcat("File Name: ", filename));
        set(t,'Interpreter','none');

        % compare imag part
        subplot(5,1,2)
        plot(imag(samples_antenna_tx_matlab_resampled(:,i)));
        hold on
        plot(imag(samples_antenna_tx_cpp(:,i)), 'r');
        legend('Matlab', 'C++');
        grid on

        % compare real before and after resampling
        subplot(5,1,3)
        plot(t_dect, real(samples_antenna_tx_matlab(:,i)));
        hold on
        plot(t_hw, real(samples_antenna_tx_cpp(:,i)), 'r');
        legend('Matlab before resampling', 'C++ resampled');
        grid on

        % compare imag before and after resampling
        subplot(5,1,4)
        plot(t_dect, imag(samples_antenna_tx_matlab(:,i)));
        hold on
        plot(t_hw, imag(samples_antenna_tx_cpp(:,i)), 'r');
        legend('Matlab before resampling', 'C++ resampled');
        grid on

        % observe frequency shift
        subplot(5,1,5)
        spec_tmp = sum(samples_antenna_tx_matlab_resampled, 2);
        spec_tmp = fftshift(fft(spec_tmp));
        spec_tmp = abs(spec_tmp).^2;
        spec_tmp = pow2db(spec_tmp);
        plot(spec_tmp);
        grid on
    end
end