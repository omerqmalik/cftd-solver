function S_obj = userplot_makeSpatialFieldMovie(cav_dir,num,beta,t0,t1,numframes)
    S_EcoeffsAll = structdata_load(cav_dir,num,'E','coeffs',1,beta);
    t = S_EcoeffsAll.t;
    
    [~,t0_ind] = helpers_getClosestMatch(t0,t);
    [~,t1_ind] = helpers_getClosestMatch(t1,t);
    
    d = fix((t1_ind - t0_ind)/(numframes - 1));
    time_inds = t0_ind:d:t1_ind;
    if time_inds(end) < t1_ind
        time_inds = [time_inds (time_inds(end)+d)];
    end
    
    S_dataArray = structdata_getDataArray(S_EcoeffsAll,'time',time_inds);
    S_obj       = plotting_plotData(S_dataArray,'SPCFLDmov','2Dcoeffs','');
end