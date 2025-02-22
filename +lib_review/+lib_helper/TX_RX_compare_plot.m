function [] = TX_RX_compare_plot(samples_antenna_tx_matlab_resampled, samples_antenna_rx_cpp, IQ_err_max_idx)

    % plot both signals
    fig = figure(1);

    % plot comparison
    x_margin = 50;
    y_margin = 50;

    % figure position
    fig_pos_x = x_margin;
    fig_pos_y = -y_margin;
    fig_pos = [fig_pos_x fig_pos_y];
    movegui(fig,fig_pos);
    clf()

    % compare real part
    subplot(2,1,1)
    plot(real(samples_antenna_tx_matlab_resampled));
    hold on
    plot(real(samples_antenna_rx_cpp), 'r');
    xline(IQ_err_max_idx)
    legend('Matlab', 'C++');
    grid on

    t = title("TX RX SISO comparison");
    set(t,'Interpreter','none');

    % compare imag part
    subplot(2,1,2)
    plot(imag(samples_antenna_tx_matlab_resampled));
    hold on
    plot(imag(samples_antenna_rx_cpp), 'r');
    xline(IQ_err_max_idx)
    legend('Matlab', 'C++');
    grid on
end