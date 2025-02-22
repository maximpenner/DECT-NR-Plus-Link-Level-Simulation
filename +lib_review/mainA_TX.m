clear all;
close all;
clc;

% In C++:
%
%   1) Disable phase rotation.
%   2) Set GI percentage to 100.
%   3) En- or disable cover sequence.
%   4) Enable JSON export.
%
% In Matlab:
%
%   1) En- or disable cover sequence.

warning on

% load and separate all json files
[filenames, n_files] = lib_util.get_all_filenames('results');
[tx_filenames, ~, ~] = lib_review.lib_helper.json_separate(filenames, n_files);

plot_every_packet = false;

% process each file
for i=1:1:numel(tx_filenames)

    % extract file folder and name
    file_struct = tx_filenames(i);
    filefolder = file_struct.folder;
    filename = file_struct.name;
    ffn = fullfile(filefolder, filename);
    disp(file_struct.name);

    % load json
    tx_json_struct = lib_review.lib_helper.json_load(ffn);

    % compare TX signals of C++ and Matlab numerically
    [~, ...
     samples_antenna_tx_matlab, ...
     samples_antenna_tx_matlab_resampled, ...
     samples_antenna_tx_cpp, ...
     dect_samp_rate, ...
     hw_samp_rate] = lib_review.lib_helper.TX_compare_numerically(ffn, tx_json_struct);

    % skip if numeric comparison failed gracefully
    if isempty(samples_antenna_tx_matlab)
        continue;
    end

    % numerical comparison successful, now plot for visual confirmation
    if plot_every_packet == true
        lib_review.lib_helper.TX_compare_plot(  filename, ...
                                                samples_antenna_tx_matlab, ...
                                                samples_antenna_tx_matlab_resampled, ...
                                                samples_antenna_tx_cpp, ...
                                                dect_samp_rate, ...
                                                hw_samp_rate);
    end
end

