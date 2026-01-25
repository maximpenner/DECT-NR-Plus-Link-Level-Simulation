function [filenames, n_files] = get_all_filenames(folderpath)

    filenames = dir(folderpath);
    
    % remove . and ..
    filenames(1) = [];
    filenames(1) = [];
    
    n_files = numel(filenames);
end
