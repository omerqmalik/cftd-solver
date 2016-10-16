function userdata_saveY(Yint,data_dir,data_id,data_type,pumpstep)
    filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep);
    save(filename,'Yint');
end