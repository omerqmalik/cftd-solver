function rawdata_save(T,Y,data_dir,data_id,data_type,pumpstep)
    filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep);
    save(filename,'T','Y');
end