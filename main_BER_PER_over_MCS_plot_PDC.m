clear all;
close all;

% This scripts plots the results generated with main_BER_PER_over_MCS.m.

load('results/var_all.mat');

ber_global = n_bits_PDC_error_global./n_bits_PDC_sent_global;
per_global = n_packets_PDC_error_global./n_packets_PDC_sent_global;

% required plot configuration
K = db2pow(9.0);

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% plot configuration
colors_vec = [0,      0.4470, 0.7410;...
              0.8500, 0.3250, 0.0980;...
              0.9290, 0.6940, 0.1250;...
              0.4940, 0.1840, 0.5560;...
              0.4660, 0.6740, 0.1880;...
              0.3010, 0.7450, 0.9330;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840;...
              0.6350, 0.0780, 0.1840];
legend_font_size = 8;
marker_size = 4;
axis_lim = [-15 55 1e-7 1e1];
legend_location = 'NorthEast';

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % extract
    mcs_index = mcs_index_vec(cnt);
    snr_db_vec = snr_db_vec_global(cnt,:);
    bps = bps_global(cnt);
    tbs = tbs_global(cnt);
    colors = colors_vec(cnt, :);
    ber = ber_global(cnt,:);

    % plot BER references curves
    lib_util.lib_plot_print.plot_ber_over_snr_references(mcs_index, snr_db_vec, bps, tbs, colors, K);

    % plot measured BER
    str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs));
    lineH = semilogy(snr_db_vec, ber_global(cnt,:), '-.o','DisplayName',str, 'Color', colors, 'MarkerSize', marker_size, 'MarkerFaceColor', colors);
    hold on
end

% configure plot
title('BER uncoded')
xlabel('SNR (dB)')
ylabel('BER uncoded')
legend('Location',legend_location, 'FontSize', legend_font_size)
grid on
axis(axis_lim)
set(gca, 'ColorOrder', jet(100))

savefig('results/A_BER_SNR_PCC.fig')

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % extract
    mcs_index = mcs_index_vec(cnt);
    snr_db_vec = snr_db_vec_global(cnt,:);
    bps = bps_global(cnt);
    tbs = tbs_global(cnt);
    colors = colors_vec(cnt, :);
    ber = ber_global(cnt,:);

    % plot BER references curves
    lib_util.lib_plot_print.plot_ber_over_ebn0_references(mcs_index, snr_db_vec, bps, tbs, colors, K);

    % convert subcarrier snr to ebn0
    EbN0_vec = db2pow(snr_db_vec)/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

    % plot measured BER
    str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs));
    lineH = semilogy(EbN0_dB_vec, ber, '-.o','DisplayName', str, 'Color', colors, 'MarkerSize', marker_size, 'MarkerFaceColor', colors);
end

% configure plot
title('BER uncoded')
xlabel('EbN0 (dB)')
ylabel('BER uncoded')
legend('Location',legend_location, 'FontSize', legend_font_size)
grid on
axis(axis_lim)
set(gca, 'ColorOrder', jet(100))

savefig('results/B_BER_EBN0_PCC.fig')

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % extract
    mcs_index = mcs_index_vec(cnt);
    colors = colors_vec(cnt, :);
    per = per_global(cnt,:);

    % plot measured PER
    str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs));
    semilogy(snr_db_vec, per,'-o','DisplayName',str, 'Color', colors, 'MarkerSize', marker_size, 'MarkerFaceColor', colors);
    hold on
end

% configure plot
title('PER')
xlabel('SNR (dB)')
ylabel('PER')
legend('Location',legend_location, 'FontSize', legend_font_size)
grid on
axis(axis_lim)
set(gca, 'ColorOrder', jet(100))

savefig('results/C_PER_SNR_PCC.fig')