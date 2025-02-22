function [tx_to_sync_points_optimal_idx] = TX_compare_sync_points(tx_filenames, string_sync_points_src, sync_fine_peak_time_64_vec, sync_point_correct_limit)

    % 0 means no optimal file found
    tx_to_sync_points_optimal_idx = ones(size(tx_filenames)) * 0;

    % compare tx times against job queue
    sync_point_err_cnt = 0;
    for i=1:1:numel(tx_filenames)
    
        % extract file folder and name
        file_struct = tx_filenames(i);
        filefolder = file_struct.folder;
        filename = file_struct.name;
        ffn = fullfile(filefolder, filename);
    
        % load json
        tx_json_struct = lib_review.lib_helper.json_load(ffn);
    
        % find the closest RX time
        [tx_to_sync_point_difference, tx_to_sync_points_optimal_idx_single] = ...
            min(abs(sync_fine_peak_time_64_vec-tx_json_struct.meta.transmission_time_64));
    
        % ignore TX file is no reasonable RX file found
        if tx_to_sync_point_difference > sync_point_correct_limit
            sync_point_err_cnt = sync_point_err_cnt + 1;
        else
            tx_to_sync_points_optimal_idx(i) = tx_to_sync_points_optimal_idx_single;
        end
    end
    
    fprintf("%s  TX packets=%d  joq_queue syncs=%d  misses=%d  SER=%f\n",   string_sync_points_src,...
                                                                            numel(tx_filenames), ...
                                                                            numel(sync_fine_peak_time_64_vec), ...
                                                                            sync_point_err_cnt, ...
                                                                            sync_point_err_cnt/numel(tx_filenames));

end