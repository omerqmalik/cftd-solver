function [T,Y] = rawdata_load(data_dir,data_id,data_type,pumpstep)
    filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep);
    load(filename,'T','Y');
end