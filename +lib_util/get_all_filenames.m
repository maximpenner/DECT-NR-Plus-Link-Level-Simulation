function [filenames, n_files] = get_all_filenames(folderpath)

    filenames = dir(folderpath);
    
    % remove . and ..
    filenames(1) = [];
    filenames(1) = [];
    
    n_files = numel(filenames);
    % if n_files == 0
    %     error("No files found!");
    % else
    %     fprintf("Read %d measurement files.\n", n_files);
    % end    
end
