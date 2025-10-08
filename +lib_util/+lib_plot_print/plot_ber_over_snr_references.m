function [] = plot_ber_over_snr_references(mcs_index, snr_db_vec, bps, tbs, colors, K)

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
    EbN0_vec = db2pow(snr_db_vec)/bps;
    EbN0_dB_vec = pow2db(EbN0_vec);

    % AWGN, diversity order 1
    if M==2
        ber_rayleigh = berawgn(EbN0_dB_vec,'psk',M,1);
    else
        ber_rayleigh = berawgn(EbN0_dB_vec,'qam',M,1);
    end
    str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs), ', AWGN');
    semilogy(snr_db_vec, ber_rayleigh,'-','DisplayName',str, 'Color', colors);
    hold on

        % RAYLEIGH, diversity order 1
        if M==2
            ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1);
        else
            ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1);
        end
        str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs), ', rayleigh');
        semilogy(snr_db_vec, ber_rayleigh,'--','DisplayName',str, 'Color', colors);
        hold on

            % RAYLEIGH, diversity order 2
            if M==2
                ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2);
            else
                ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2);
            end
            str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs), ', rayleigh div=2');
            semilogy(snr_db_vec, ber_rayleigh,'-.','DisplayName',str, 'Color', colors);
            hold on

                % RICIAN, diversity order 1
                if M==2
                    ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,1,K);
                else
                    ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,1,K);
                end
                str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs), ', rician');
                semilogy(snr_db_vec, ber_rayleigh,'-o','DisplayName',str, 'Color', colors);
                hold on
    
                    % RICIAN, diversity order 1
                    if M==2
                        ber_rayleigh = berfading(EbN0_dB_vec,'psk',M,2,K);
                    else
                        ber_rayleigh = berfading(EbN0_dB_vec,'qam',M,2,K);
                    end
                    str = append('MCS=',num2str(mcs_index),', TBS=',num2str(tbs), ', rician div=2');
                    semilogy(snr_db_vec, ber_rayleigh,'-d','DisplayName',str, 'Color', colors);
                    hold on
end