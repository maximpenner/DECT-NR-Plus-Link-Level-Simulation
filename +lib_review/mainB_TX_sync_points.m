clear all;
close all;
clc;

% script to compare TX and sync points processed in job_queue

warning on

% how many samples deviation are acceptable?
sync_point_correct_limit = 0;

% load and separate all json files
[filenames, n_files] = lib_util.get_all_filenames('results');
[tx_filenames, joq_queue_filename, rx_synced_filenames] = lib_review.lib_helper.json_separate(filenames, n_files);

%% sync points from joq queue
file_struct = joq_queue_filename;
filefolder = file_struct.folder;
filename = file_struct.name;
ffn = fullfile(filefolder, filename);
joq_queue_json_struct = lib_review.lib_helper.json_load(ffn);
joq_queue_sync_fine_peak_time_64_vec = joq_queue_json_struct.sync_fine_peak_time_64_vec;

tx_to_sync_points_optimal_idx_A = lib_review.lib_helper.TX_compare_sync_points( tx_filenames,...
                                                                                "Sync points from job queue:", ...
                                                                                joq_queue_sync_fine_peak_time_64_vec, ...
                                                                                sync_point_correct_limit);

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

tx_to_sync_points_optimal_idx_B = lib_review.lib_helper.TX_compare_sync_points( tx_filenames, ...
                                                                                "Sync points from rx files: ", ...
                                                                                rx_synced_sync_fine_peak_time_64_vec, ...
                                                                                sync_point_correct_limit);
