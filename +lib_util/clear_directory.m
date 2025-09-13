function [] = clear_directory(directory_path)

    % first find all recorded files
    [filenames, n_files] = lib_util.get_all_filenames(directory_path);

    % delete all files
    for i=1:1:n_files
        full_filepath = fullfile(filenames(i).folder,filenames(i).name);
        delete(full_filepath);
    end
end
