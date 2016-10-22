function rawdata_cleanAll(data_dir,data_id,data_type)
    [cT,cY,pstep] = rawdata_loadAll(data_dir,data_id,data_type);
    t             = rawdata_interpT(cT);
    
    for i = 1:length(cY)
        Y = cY{i};
        T = cT{i};
        Yint = rawdata_interpY(T,Y,t);
        userdata_saveY(Yint,data_dir,data_id,data_type,pstep(i));
    end
    userdata_saveT(t,data_dir,data_id,data_type);
end