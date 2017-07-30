function S_structdata = structdata_loadMicro(cav_dir,data_dir,num,data_id,data_type,calc_times,pstep,x0)
    S_setupdata = setup_loadParameters(cav_dir,num);
    [t,Y]  = userdata_load(data_dir,data_id,data_type,pstep);

    S_structdata   = structdata_getGeneralProperties(data_id,data_type,S_setupdata,calc_times,pstep);
    S_structdata.t = t;
    S_structdata.Y = Y;
    S_structdata.psteps = pstep;
    S_structdata.x0 = x0;
end