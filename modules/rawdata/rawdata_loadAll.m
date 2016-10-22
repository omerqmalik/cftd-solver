function [cT,cY,fnum] = rawdata_loadAll(data_dir,data_id,data_type)
    datafiles = dir([data_dir '/' data_id data_type '_*.mat']);
    [datafiles,fnum] = rawdata_orderFilesNumerically(datafiles);
    len_DF = length(datafiles);
    
    cT = cell(1,len_DF);
    cY = cell(1,len_DF);
    for i = 1:len_DF
        load([data_dir '/' datafiles(i).name]);

        cT{i} = T;
        cY{i} = Y;
    end
end