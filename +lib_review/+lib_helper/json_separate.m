function [tx_filenames, joq_queue_filename, rx_synced_filenames] = json_separate(filenames, n_files)

    % sanity check
    if numel(filenames) ~= n_files
        error("Length of struct not n_files.");
    end

    tx_filenames = filenames;
    joq_queue_filename = filenames;
    rx_synced_filenames = filenames;

    % go over each file, starting with the last one
    for i = n_files : -1 : 1

        % remove entry in tx_filenames if not starting with the correct letters
        if strncmpi(filenames(i).name, 'tx_packet', numel('tx_packet')) == true
            joq_queue_filename(i) = [];
            rx_synced_filenames(i) = [];
        elseif strncmpi(filenames(i).name, 'job_queue', numel('job_queue')) == true
            tx_filenames(i) = [];
            rx_synced_filenames(i) = [];
        elseif strncmpi(filenames(i).name, 'rx_packet', numel('rx_packet')) == true
            tx_filenames(i) = [];
            joq_queue_filename(i) = [];
        else
            tx_filenames(i) = [];
            joq_queue_filename(i) = [];
            rx_synced_filenames(i) = [];
        end
    end
end