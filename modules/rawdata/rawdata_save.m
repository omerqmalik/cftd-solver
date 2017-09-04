function rawdata_save(T,Y,data_dir,data_id,data_type,pumpstep,opt_checkpoint)
    if nargin == 6
        filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep);
    elseif nargin == 7
        filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep,opt_checkpoint);
    end
    save(filename,'T','Y');
end