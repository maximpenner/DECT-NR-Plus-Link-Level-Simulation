function json_struct = json_load(folder_filename)

    fid = fopen(folder_filename);
    raw = fread(fid,inf);
    fclose(fid);

    % interpret
    str = char(raw');
    json_struct = jsondecode(str);
end