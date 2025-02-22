function [] = json_save(folder_filename, json_str)

    fid = fopen(folder_filename,'wt');
    fprintf(fid, json_str);
    fclose(fid);
    
end