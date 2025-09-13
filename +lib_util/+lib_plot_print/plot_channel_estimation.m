function [] = plot_channel_estimation(ch_estim, n_guards_bottom, n_guards_top)

    N_RX = numel(ch_estim);

    ch_estim_single = ch_estim{1};

    [N_subc, N_symb, N_TS] =  size(ch_estim_single);

    subc_idx = -N_subc/2 : 1 : (N_subc/2-1);

    figure()
    clf()

    for i=1:N_RX

        % extract
        ch_estim_single = ch_estim{i};

        for j=1:N_TS
        
            % extract
            ch_estim_single_TS = ch_estim_single(:,:,j);

            % prepare subplot
            subplot(N_RX, N_TS, (i-1)*N_TS + j);
            hold on;

            % fill subplot
            for k=1:N_symb

                % extract
                ch_estim_single_TS_symb = ch_estim_single_TS(:, k);

                ch_estim_single_TS_symb(1:n_guards_bottom) = 0;
                ch_estim_single_TS_symb(end-n_guards_top+1:end) = 0;

                curve2plot = abs(ch_estim_single_TS_symb);
                curve2plot = mag2db(curve2plot);

                % plot
                plot(subc_idx, curve2plot);
            end

            % finalize subplot
            title_str = "Ch Estim: RX " + num2str(i-1) + " TS " + num2str(j-1) + " for " + num2str(N_symb) + " Symbols";
            title(title_str);
            grid on
            grid minor
            xlim([-N_subc/2 N_subc/2])
            ylim([-30 10])
            xlabel("Subcarrier Index")
            ylabel("Absolute dB")
        end
    end
end
