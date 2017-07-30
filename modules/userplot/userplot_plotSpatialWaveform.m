function S_obj = userplot_plotSpatialWaveform(cav_dir,num,x0,beta,t0,t1,numframes)
    if nargin < 5
        t1 = t0;
        numframes = 1;
    end
    S_EcoeffsAll = structdata_load(cav_dir,num,'E','coeffs',1,x0,beta); 
    
    time_inds = helpers_getSptlFldMovTimeIndices(S_EcoeffsAll,t0,t1,numframes);
    
    S_EcoeffsAll.t_ind = time_inds;
    S_EcoeffsAll.Y = S_EcoeffsAll.Y(time_inds,:);
        
    S_obj = plotting_plotData(S_EcoeffsAll,'SPCFLD','2Dcoeffs','');
end