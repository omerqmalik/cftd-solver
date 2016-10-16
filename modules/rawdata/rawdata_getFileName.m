function filename = rawdata_getFileName(data_dir,data_id,data_type,num)
    if strcmp(num,'t')
        filename = [data_dir,'/',data_id,data_type,'_t','.mat'];
    else
        filename = [data_dir,'/',data_id,data_type,'_',num2str(num),'.mat'];
    end
end