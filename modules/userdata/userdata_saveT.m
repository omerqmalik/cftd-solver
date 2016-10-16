function userdata_saveT(t,data_dir,data_id,data_type)
    filename = rawdata_getFileName(data_dir,data_id,data_type,'t');
    save(filename,'t');
end