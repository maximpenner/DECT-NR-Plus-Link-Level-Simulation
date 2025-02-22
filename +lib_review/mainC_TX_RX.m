clear all;
close all;
clc;

% script to compare TX and RX of C++ (exported with loopback firmware) with Matlab encoding/decoding

warning on

% how many samples deviation are acceptable?
sync_point_correct_limit = 2;

% load and separate all json files
[filenames, n_files] = lib_util.get_all_filenames('results');
[tx_filenames, ~, rx_synced_filenames] = lib_review.lib_helper.json_separate(filenames, n_files);

%% sync points from rx files
rx_synced_sync_fine_peak_time_64_vec = zeros(numel(rx_synced_filenames), 1);
for i=1:1:numel(rx_synced_filenames)

    % extract file folder and name
    file_struct = rx_synced_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);

    % load json
    rx_json_struct = lib_review.lib_helper.json_load(ffn);

    rx_synced_sync_fine_peak_time_64_vec(i) = rx_json_struct.meta.sync_fine_peak_time_64;
end

% find mapping between tx times and sync points
tx_to_sync_points_optimal_idx = lib_review.lib_helper.TX_compare_sync_points(   tx_filenames, ...
                                                                                "Sync points from rx files: ", ...
                                                                                rx_synced_sync_fine_peak_time_64_vec, ...
                                                                                sync_point_correct_limit);

%% check every tx packet numerically, and then compare with rx packet

results = cell(numel(tx_filenames), 1);

plot_every_packet = false;

for i=1:1:numel(tx_filenames)

    if tx_to_sync_points_optimal_idx(i) <= 0
        continue;
    end

    % extract file folder and name
    file_struct = tx_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);
    disp(file_struct.name);

    % load tx json
    tx_json_struct = lib_review.lib_helper.json_load(ffn);

    % compare TX signals of C++ and Matlab numerically
    [tx, ...
     samples_antenna_tx_matlab, ...
     samples_antenna_tx_matlab_resampled, ...
     samples_antenna_tx_cpp, ...
     dect_samp_rate, ...
     hw_samp_rate] = lib_review.lib_helper.TX_compare_numerically(ffn, tx_json_struct);

    % extract RX file folder and name
    file_struct = rx_synced_filenames(tx_to_sync_points_optimal_idx(i));
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);

    % load rx json
    rx_json_struct = lib_review.lib_helper.json_load(ffn);

    % first compare packets numerically
    result_signal = lib_review.lib_helper.TX_RX_compare_numerically(tx, samples_antenna_tx_matlab_resampled, tx_json_struct, rx_json_struct);

    % save
    results(i) = {result_signal};

    % then plot
    if plot_every_packet == true
        lib_review.lib_helper.TX_RX_compare_plot(samples_antenna_tx_matlab_resampled, result_signal.samples_antenna_rx_cpp, result_signal.IQ_err_max_idx);
    end
end

figure(1)
clf()
plot_snr_hist(results, "IQ_snr_dB");

figure(2)
clf()
subplot(2,1,1)
plot_snr_hist(results, "PCC_snr_dB");
subplot(2,1,2)
plot_ber_hist(results, "BER_PCC")

figure(3)
clf()
subplot(2,1,1)
plot_snr_hist(results, "PDC_snr_dB");
subplot(2,1,2)
plot_ber_hist(results, "BER_PDC")

function plot_snr_hist(results, fieldname)
    hist_container = get_histogramm_of_field(results, fieldname);
    histogram(hist_container);
    xlim([0 max(15, 1.5*max(hist_container))])
    title_str = fieldname + " for " + num2str(numel(results)) +  " files";
    title(title_str, 'Interpreter', 'none');
end

function plot_ber_hist(results, fieldname)
    hist_container = get_histogramm_of_field(results, fieldname);
    histogram(hist_container, 200);
    xlim([-0.2 1.2])
    title_str = fieldname + " for " + num2str(numel(results)) +  " files";
    title(title_str, 'Interpreter', 'none');
end

function hist_container = get_histogramm_of_field(results, fieldname)
    hist_container = zeros(numel(results), 1);
    for i = numel(results):-1:1
        if isempty(results(i)) == false
            S = results{i};
            hist_container(i) = extractfield(S, fieldname);
        else
            hist_container(i) = [];
        end
    end
end
