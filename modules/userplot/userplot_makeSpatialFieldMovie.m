function S_obj = userplot_makeSpatialFieldMovie(cav_dir,num,beta,t0,t1,numframes)
    S_EcoeffsAll = structdata_load(cav_dir,num,'E','coeffs',1,beta);
    time_inds = helpers_getSptlFldMovTimeIndices(S_EcoeffsAll,t0,t1,numframes);
    
    S_dataArray = structdata_getDataArray(S_EcoeffsAll,'time',time_inds);
    S_obj       = plotting_plotData(S_dataArray,'SPCFLDmov','2Dcoeffs','');
end