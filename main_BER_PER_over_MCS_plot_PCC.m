clear all;
close all;

% This scripts plots the results generated with main_BER_PER_over_MCS.m.

load('results/var_all.mat');

ber_global = n_bits_PCC_error_global./n_bits_PCC_sent_global;
per_global = n_packets_PCC_error_global./n_packets_PCC_sent_global;

% overwrite, fixed for PCC, always QPSK
bps_global = ones(size(bps_global))*2;

if tx.tx_config.PLCF_type == 1
    plcf_length = 40;
elseif tx.tx_config.PLCF_type == 2
    plcf_length = 80;
else
    error('Unknown PLCF length.');
end

tbs_global = ones(size(tbs_global))*plcf_length;

% required plot configuration
K = db2pow(9.0);

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

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

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % bits per symbol
    bps = bps_global(cnt);
    M = 2^bps;

    % convert subcarrier snr to ebn0
    %
    %   S/N = R_b * E_b / (B * N0)
    %
    %   R_b = bps / T
    %   B = 1 / T
    %
    %   S/N = bps * E_b/N0
    %
    % source: https://www.gaussianwaves.com/2008/11/relation-between-ebn0-and-snr-2/
    EbN0_vec = db2pow(snr_db_vec_global(cnt,:))/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

    % AWGN, diversity order 1
    if M==2
        ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,1);
    else
        ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M,1);
    end
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', AWGN');
    semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-','DisplayName',str, 'Color', colors(cnt,:));
    hold on

        % RAYLEIGH, diversity order 1
        if M==2
            ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
        else
            ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
        end
        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh');
        semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'--','DisplayName',str, 'Color', colors(cnt,:));
        hold on

            % RAYLEIGH, diversity order 2
            if M==2
                ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
            else
                ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
            end
            str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh div=2');
            semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-.','DisplayName',str, 'Color', colors(cnt,:));
            hold on

                % RICIAN, diversity order 1
                if M==2
                    ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
                else
                    ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
                end
                str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician');
                semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-o','DisplayName',str, 'Color', colors(cnt,:));
                hold on
    
                    % RICIAN, diversity order 1
                    if M==2
                        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
                    else
                        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
                    end
                    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician div=2');
                    semilogy(snr_db_vec_global(cnt,:), ber_rayleigh,'-d','DisplayName',str, 'Color', colors(cnt,:));
                    hold on

    % SIMULATION

    % our measured ber
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    lineH = semilogy(snr_db_vec_global(cnt,:), ber_global(cnt,:), '-.o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on

    % general settings
    title('BER uncoded')
    xlabel('SNR (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
end

savefig('results/A_BER_SNR_PCC.fig')

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)

    % bits per symbol
    bps = bps_global(cnt);
    M = 2^bps;

    % convert subcarrier snr to ebn0
    %
    %   S/N = R_b * E_b / (B * N0)
    %
    %   R_b = bps / T
    %   B = 1 / T
    %
    %   S/N = bps * E_b/N0
    %
    % source: https://www.gaussianwaves.com/2008/11/relation-between-ebn0-and-snr-2/
    EbN0_vec = db2pow(snr_db_vec_global(cnt,:))/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

    % AWGN, diversity order 1
    if M==2
        ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,1);
    else
        ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M,1);
    end
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', AWGN');
    semilogy(EbN0_dB_vec, ber_rayleigh,'-','DisplayName',str, 'Color', colors(cnt,:));
    hold on

        % RAYLEIGH, diversity order 1
        if M==2
            ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
        else
            ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
        end
        str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh');
        semilogy(EbN0_dB_vec, ber_rayleigh,'--','DisplayName',str, 'Color', colors(cnt,:));
        hold on

            % RAYLEIGH, diversity order 2
            if M==2
                ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
            else
                ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
            end
            str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rayleigh div=2');
            semilogy(EbN0_dB_vec, ber_rayleigh,'-.','DisplayName',str, 'Color', colors(cnt,:));
            hold on

                % RICIAN, diversity order 1
                if M==2
                    ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
                else
                    ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
                end
                str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician');
                semilogy(EbN0_dB_vec, ber_rayleigh,'-o','DisplayName',str, 'Color', colors(cnt,:));
                hold on
    
                    % RICIAN, diversity order 1
                    if M==2
                        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
                    else
                        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
                    end
                    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)), ', rician div=2');
                    semilogy(EbN0_dB_vec, ber_rayleigh,'-d','DisplayName',str, 'Color', colors(cnt,:));
                    hold on

    % SIMULATION

    % our measured ber
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    lineH = semilogy(EbN0_dB_vec, ber_global(cnt,:), '-.o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on

    % general settings
    title('BER uncoded')
    xlabel('EbN0 (dB)')
    ylabel('BER uncoded')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
    set(gca, 'ColorOrder', jet(100))
end

savefig('results/B_BER_EBN0_PCC.fig')

figure()
clf()
for cnt = 1:1:numel(mcs_index_vec)
    str = append('MCS=',num2str(mcs_index_vec(cnt)),', TBS=',num2str(tbs_global(cnt)));
    semilogy(snr_db_vec_global(cnt,:), per_global(cnt,:),'-o','DisplayName',str, 'Color', colors(cnt,:), 'MarkerSize', marker_size, 'MarkerFaceColor', colors(cnt,:));
    hold on
    title('Packet Decoding Performance')
    xlabel('SNR (dB)')
    ylabel('PER')
    legend('Location',legend_location, 'FontSize', legend_font_size)
    grid on
    axis(axis_lim)
end

savefig('results/C_PER_SNR_PCC.fig')
