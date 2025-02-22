clear all;
close all;
clc;

% script to plot PER results from C++

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

    nof_packets = json_struct.nof_packets;
    MCS = json_struct.MCS_index;
    TB_bits= json_struct.TB_bits;

    snr_vec = json_struct.data.snr_vec;
    PER_pcc_crc= json_struct.data.PER_pcc_crc;
    PER_pcc_crc_content= json_struct.data.PER_pcc_crc_content;
    PER_pdc_crc = json_struct.data.PER_pdc_crc;

    % plot FEC curves into specific figure

    x_margin_px = 50;
    y_margin_px = 150;
    fig_corner_spacing_in_px = 25;

    fig_nr = 1 + MCS;
    fig_pos_x = x_margin_px + MCS * fig_corner_spacing_in_px;
    fig_pos_y = -(y_margin_px + MCS * fig_corner_spacing_in_px);

    % create and move figure window
    fig = figure(fig_nr);
    movegui(fig, [fig_pos_x fig_pos_y]);

    title_str = "MCS=" + num2str(MCS) + " TB_bits=" + TB_bits;

    % plot SNR comparison, should be a straight line
    semilogy(snr_vec, PER_pcc_crc, '-*', 'DisplayName', "PER_pcc_crc");
    hold on
    semilogy(snr_vec, PER_pcc_crc_content, '-*', 'DisplayName', "PER_pcc_crc_content");
    semilogy(snr_vec, PER_pdc_crc, '-*', 'DisplayName', "PER_pdc_crc");
    grid on
    hold on
    xlim([-5 20]);
    ylim([1e-6 5*1e2]);
    xlabel('SNR')
    ylabel('SNR Measured')

    t = title(title_str);
    set(t,'Interpreter','none');
    l = legend('Location', 'southeast');
    set(l,'Interpreter','none');

    pause(0.5);
end