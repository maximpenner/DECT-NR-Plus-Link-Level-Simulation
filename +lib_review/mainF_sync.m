clear all;
%close all;
clc;

% script to plot PER results from C++

warning on

% load all json filenames
[filenames, n_files] = lib_util.get_all_filenames('results');

figure(1)
clf();

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

    nof_antennas = json_struct.nof_antennas;
  
    waveform_power = json_struct.waveform_power;
    waveform_metric = json_struct.waveform_metric;
    waveform_metric_smooth = json_struct.waveform_metric_smooth;
    waveform_metric_max_idx = json_struct.waveform_metric_max_idx;
    metric_smoother_length = json_struct.metric_smoother_length;
    metric_smoother_bos_offset_to_center_samples = json_struct.metric_smoother_bos_offset_to_center_samples;

    % revert concatenation
    waveform_power = reshape(waveform_power, [], nof_antennas);
    waveform_metric = reshape(waveform_metric, [], nof_antennas);
    waveform_metric_smooth = reshape(waveform_metric_smooth, [], nof_antennas);

    figure(1)
    clf();

    n_elem_per_antenna = size(waveform_power,1);

    time_axis = 1:1:n_elem_per_antenna;

    for i=1:nof_antennas

        % C++ starts indexing at 0
        metric_idx_max_matlab_cpp_equivalent = waveform_metric_max_idx(i) + 1;

        % find index of local maximum
        [~, metric_idx_max_matlab] = max(waveform_metric(:,i));
        [~, metric_smooth_idx_max_matlab] = max(waveform_metric_smooth(:,i));

        %assert(metric_idx_max_matlab_cpp_equivalent == metric_idx_max_matlab);

        % plot power
        subplot(nof_antennas, 2, (i-1)*2 + 1);
        plot(time_axis, waveform_power(:,i))
        hold on
        xline(metric_idx_max_matlab_cpp_equivalent);
        xline(metric_smooth_idx_max_matlab, 'r');
        grid on
        xlim([-5 n_elem_per_antenna+5]);
        %ylim([-0.01 0.01]);

        % plot metric
        subplot(nof_antennas, 2, (i-1)*2 + 2);
        plot(time_axis, waveform_metric(:,i))
        hold on
        plot(time_axis, waveform_metric_smooth(:,i), 'k')
        xline(metric_idx_max_matlab_cpp_equivalent);
        xline(metric_smooth_idx_max_matlab, 'r');
        xline(metric_idx_max_matlab, 'g');
        grid on
        xlim([-5 n_elem_per_antenna+5]);
        ylim([-0.5 1.5]);
    end

    1;
end