function userdata_saveY(Yint,data_dir,data_id,data_type,pstep)
    filename = rawdata_getFileName(data_dir,data_id,data_type,pstep);
    save(filename,'Yint');
end