function filename = rawdata_getFileName(data_dir,data_id,data_type,pstep,opt_checkpoint)
    if strcmp(pstep,'t')
        filename = [data_dir,'/',data_id,data_type,'_t','.mat'];
    else
        if nargin == 4
            filename = [data_dir,'/',data_id,data_type,'_',num2str(pstep),'.mat'];
        elseif nargin == 5
            filename = [data_dir,'/',data_id,data_type,'_',num2str(pstep),'_',num2str(opt_checkpoint),'.mat'];
        end
    end
end