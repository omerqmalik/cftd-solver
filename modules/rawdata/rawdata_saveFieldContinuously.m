function rawdata_saveFieldContinuously(T,Y_field,data_dir,data_id,data_type,pstep,tstep)
    Yfield_next = Y_field;
    T_next   = T;

    if tstep == 1   %Always create a new file if first time-step (ie: overwrite existing file)
        rawdata_save(T_next,Yfield_next,data_dir,data_id,data_type,pstep);
    elseif exist(rawdata_getFileName(data_dir,data_id,data_type,pstep),'file') == 2
        [T_this,Yfield_this] = rawdata_load(data_dir,data_id,data_type,pstep);
        Yfield_this = [Yfield_this; Yfield_next];
        T_this      = [T_this; T_next];
        rawdata_save(T_this,Yfield_this,data_dir,data_id,data_type,pstep);
    else
        rawdata_save(T_next,Yfield_next,data_dir,data_id,data_type,pstep);
    end
end