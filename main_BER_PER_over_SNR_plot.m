clear all;
close all;

load('results/var_all.mat');

% PCC
ber = n_bits_PCC_error./n_bits_PCC_sent;
per = n_packets_PCC_error./n_packets_PCC_sent;
plot_BER_PER_over_MCS("PCC", mcs, snr_db, ones(size(bps))*2, ones(size(tbs))*tx.tx_config.PLCF_type*40, ber, per);

% PDC
ber = n_bits_PDC_error./n_bits_PDC_sent;
per = n_packets_PDC_error./n_packets_PDC_sent;
plot_BER_PER_over_MCS("PDC", mcs, snr_db, bps, tbs, ber, per);

function [] = plot_BER_PER_over_MCS(prefix, mcs, snr_db, bps, tbs, ber, per)

    % plot configuration
    colors = [0,      0.4470, 0.7410;...
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

    % K-factor of Rician channel
    K = db2pow(9.0);

    % BER over EbN0
    figure()
    clf()
    for cnt = 1:1:numel(mcs)
        plot_ber_references("EbN0", mcs(cnt), snr_db, bps(cnt), tbs(cnt), colors(cnt, :), K);
    
        EbN0_vec = db2pow(snr_db)/bps(cnt);
        EbN0_dB_vec = pow2db(EbN0_vec);
        
        str = append('MCS=', num2str(mcs(cnt)), ', TBS=', num2str(tbs(cnt)));
        semilogy(EbN0_dB_vec, ber(cnt,:), '-.o','DisplayName', str, 'Color', colors(cnt, :), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt, :));
    end

    title(prefix)
    xlabel('EbN0 (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
    savefig("results/" + prefix + "_BER_EBN0.fig")
    
    % BER over SNR
    figure()
    clf()
    for cnt = 1:1:numel(mcs)
        plot_ber_references("SNR", mcs(cnt), snr_db, bps(cnt), tbs(cnt), colors(cnt, :), K);
    
        str = append('MCS=', num2str(mcs(cnt)), ', TBS=', num2str(tbs(cnt)));
        semilogy(snr_db, ber(cnt,:), '-.o','DisplayName',str, 'Color', colors(cnt, :), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt, :));
        hold on
    end

    title(prefix)
    xlabel('SNR (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
    savefig("results/" + prefix + "_BER_SNR.fig")
    
    % PER
    figure()
    clf()
    for cnt = 1:1:numel(mcs)
        str = append('MCS=', num2str(mcs(cnt)), ', TBS=', num2str(tbs(cnt)));
        semilogy(snr_db, per(cnt,:),'-o','DisplayName',str, 'Color', colors(cnt, :), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt, :));
        hold on
    end

    title(prefix)
    xlabel('SNR (dB)')
    ylabel('PER')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
    savefig("results/" + prefix + "_PER_SNR.fig")
end

function [] = plot_ber_references(reference, mcs, snr_db, bps, tbs, colors, K)
    % convert subcarrier snr to EbN0
    %
    %   S/N = R_b * E_b / (B * N0)
    %
    %   R_b = bps / T
    %   B = 1 / T
    %
    %   S/N = bps * E_b/N0
    %
    % source: https://www.gaussianwaves.com/2008/11/relation-between-ebn0-and-snr-2/
    EbN0_vec = db2pow(snr_db)/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

    M = 2^bps;

    if strcmp(reference, "SNR")
        reference = snr_db;
    elseif strcmp(reference, "EbN0")
        reference = EbN0_dB_vec;
    else
        error("unknown reference");
    end

    % AWGN, diversity order 1
    if M==2
        ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,1);
    else
        ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M,1);
    end
    str = append('MCS=',num2str(mcs),', TBS=',num2str(tbs), ', AWGN');
    semilogy(reference, ber_rayleigh,'-','DisplayName',str, 'Color', colors);
    hold on

    % Rayleigh diversity order 1
    if M==2
        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
    else
        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
    end
    str = append('MCS=',num2str(mcs),', TBS=',num2str(tbs), ', rayleigh');
    semilogy(reference, ber_rayleigh,'--','DisplayName',str, 'Color', colors);

    % Rayleigh diversity order 2
    if M==2
        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
    else
        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
    end
    str = append('MCS=',num2str(mcs),', TBS=',num2str(tbs), ', rayleigh div=2');
    semilogy(reference, ber_rayleigh,'-.','DisplayName',str, 'Color', colors);

    % Rician diversity order 1
    if M==2
        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
    else
        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
    end
    str = append('MCS=',num2str(mcs),', TBS=',num2str(tbs), ', rician');
    semilogy(reference, ber_rayleigh,'-o','DisplayName',str, 'Color', colors);
    
    % Rician diversity order 1
    if M==2
        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
    else
        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
    end
    str = append('MCS=',num2str(mcs),', TBS=',num2str(tbs), ', rician div=2');
    semilogy(reference, ber_rayleigh,'-d','DisplayName',str, 'Color', colors);
end
