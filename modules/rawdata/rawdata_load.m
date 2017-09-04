function [T,Y] = rawdata_load(data_dir,data_id,data_type,pumpstep,opt_checkpoint)
    if nargin == 4
        filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep);
    elseif nargin == 5
        filename = rawdata_getFileName(data_dir,data_id,data_type,pumpstep,opt_checkpoint);
    end
    load(filename,'T','Y');
end