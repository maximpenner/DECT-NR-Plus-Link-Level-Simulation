clear all;
close all;
clc;

% script to plot Turbo Decoder results from C++

warning on

% load all json filenames
[filenames, n_files] = lib_util.get_all_filenames('results');

% process each file
for i=1:1:n_files

    % extract file folder and name
    file_struct = filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    disp(file_struct.name)

    plot_single_file(filefolder, filename)
end

function plot_single_file(filefolder, filename)

    ffn = fullfile(filefolder, filename);
    
    json_struct = lib_review.lib_helper.json_load(ffn);

    identifier                  = json_struct.identifier;
    radio_device_class_string   = json_struct.radio_device_class_string;
    PacketLength                = json_struct.PacketLength;
    Z                           = json_struct.Z;
    MCS                         = json_struct.MCS;
    N_HARQ_RETX                 = json_struct.N_HARQ_RETX;
    N_TB_bits                   = json_struct.N_TB_bits;
    N_PDC_subc                  = json_struct.N_PDC_subc;
    G                           = json_struct.G;
    N_PACKETS                   = json_struct.N_PACKETS;

    % measure data
    SNR_dB_vec                  = json_struct.data.SNR_dB_vec;
    SNR_dB_measured_vec         = json_struct.data.SNR_dB_measured_vec;
    BER_uncoded_vec             = json_struct.data.BER_uncoded_vec;
    PER_vec                     = json_struct.data.PER_vec;

    % plot FEC curves into specific figure

    x_margin_px = 50;
    y_margin_px = 150;
    fig_corner_spacing_in_px = 25;
    MCS_max = 11;

    % group figures by Z
    if Z == 2048
        fig_nr = 1 + MCS;
        fig_pos_x = x_margin_px + MCS * fig_corner_spacing_in_px;               % margin from left side
    elseif Z == 6144
        fig_nr = 101 + MCS;
        fig_pos_x = -(x_margin_px + (MCS_max-MCS) * fig_corner_spacing_in_px);  % margin from right side
    else
        error("Unknown Z.");
    end
    fig_pos_y = -(y_margin_px + MCS * fig_corner_spacing_in_px);                % margin from top

    % create and move figure window
    fig = figure(fig_nr);
    movegui(fig, [fig_pos_x fig_pos_y]);

    title_str = "Z=" + num2str(Z) + " MCS=" + num2str(MCS);
    curve_str = "PacketLength=" + num2str(PacketLength) + " N_TB_bits=" + num2str(N_TB_bits) + " N_HARQ_RETX=" + num2str(N_HARQ_RETX);

    % plot SNR comparison, should be a straight line
    subplot(4,1,1)
    plot(SNR_dB_measured_vec, SNR_dB_measured_vec, '-*', 'DisplayName', curve_str);
    grid on
    hold on
    xlim([-20 30]);
    ylim([-20 30]);
    xlabel('SNR')
    ylabel('SNR Measured')

    t = title(title_str);
    set(t,'Interpreter','none');
    l = legend('Location', 'southeast');
    set(l,'Interpreter','none');

    % plot BER uncoded
    subplot(4,1,2)
    semilogy(SNR_dB_vec, BER_uncoded_vec, '-*', 'DisplayName', curve_str);
    grid on
    hold on
    xlim([-20 30]);
    ylim([1e-6 1e1]);
    xlabel('SNR')
    ylabel('BER uncoded')

    % plot PER, two subfigures: one with *, one without
    for j=3:4
        subplot(4,1,j)
        if j==3
            semilogy(SNR_dB_vec, PER_vec, '-*', 'DisplayName', curve_str);
        else
            semilogy(SNR_dB_vec, PER_vec, 'DisplayName', curve_str);
        end
        grid on
        hold on
        xlim([-20 30]);
        ylim([1e-6 1e1]);
        xlabel('SNR')
        ylabel('PER')
    end

    % sanity check: PER should start with 1 and must end with 0
    if(PER_vec(1) ~= 1)
        fprintf("\tFirst PER not one.\n");
    end
    if(PER_vec(end) ~= 0)
        error("\tFinal PER not zero.\n");
    end

    % sanity check: PER values must be not strictly descending
    PER_vec_sorted = sort(PER_vec, 'descend');
    err = sum(abs(PER_vec - PER_vec_sorted));
    if(err > 0)
        fprintf("\tPER not descending.\n");
        fprintf("\tFileame:      %s\n", filename);
        fprintf("\tPacketLength:  %d\n", PacketLength);
        fprintf("\tZ:             %d\n", Z);
        fprintf("\tMCS:           %d\n", MCS);
        fprintf("\tN_HARQ_RETX:   %d\n", N_HARQ_RETX);
        fprintf("\tN_TB_bits:     %d\n", N_TB_bits);
        fprintf("\tPER_vec:       ");

        % print PER vector horizontally
        for j=1:1:numel(PER_vec)
            if(PER_vec(j) == 0 || PER_vec(j) == 1)
                fprintf("%.0f ", PER_vec(j));
            else
                fprintf("%.4f ", PER_vec(j));
            end
        end
        fprintf("\n");
    end
end