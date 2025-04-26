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

    nof_antennas_limited = json_struct.nof_antennas_limited;
  
    waveform_power = json_struct.waveform_power;
    waveform_rms = json_struct.waveform_rms;
    waveform_metric = json_struct.waveform_metric;
    waveform_metric_smooth = json_struct.waveform_metric_smooth;
    waveform_metric_max_idx = json_struct.waveform_metric_max_idx;
    metric_smoother_length = json_struct.metric_smoother_length;
    metric_smoother_bos_offset_to_center_samples = json_struct.metric_smoother_bos_offset_to_center_samples;

    % revert concatenation
    waveform_power = reshape(waveform_power, [], nof_antennas_limited);
    waveform_metric = reshape(waveform_metric, [], nof_antennas_limited);
    waveform_metric_smooth = reshape(waveform_metric_smooth, [], nof_antennas_limited);

    figure(1)
    clf();

    n_elem_per_antenna = size(waveform_power,1);

    time_axis = 1:1:n_elem_per_antenna;

    for i=1:nof_antennas_limited

        % C++ starts indexing at 0
        waveform_metric_max_idx_matlab_equivalent = waveform_metric_max_idx(i) + 1;

        % plot power
        subplot(nof_antennas_limited, 2, (i-1)*2 + 1);
        plot(time_axis, waveform_rms(:,i))
        hold on
        xline(waveform_metric_max_idx_matlab_equivalent);
        grid on
        xlim([-5 n_elem_per_antenna+5]);

        % plot metric
        subplot(nof_antennas_limited, 2, (i-1)*2 + 2);
        plot(time_axis, waveform_metric(:,i))
        hold on
        xline(waveform_metric_max_idx_matlab_equivalent);
        grid on
        xlim([-5 n_elem_per_antenna+5]);
        ylim([-0.5 1.5]);

        % plot smoothed metric
        plot(time_axis, waveform_metric_smooth(:,i), 'k')
    end

    1;
end