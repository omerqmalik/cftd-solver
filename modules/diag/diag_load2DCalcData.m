function S_diagdata = diag_load2DCalcData(cav_dir,num,data_id,data_type,pstep)
    S_setupdata = setup_loadParameters(cav_dir,num);
    calc_times  = benchmark_loadTimes([cav_dir '/' S_setupdata.times_dir]);
    
    [t,Y]  = rawdata_load([cav_dir '/' S_setupdata.data_dir],data_id,data_type,pstep);

    S_diagdata   = structdata_getGeneralProperties(data_id,data_type,S_setupdata,calc_times,pstep);
    S_diagdata.t = t;
    S_diagdata.Y = Y;
end