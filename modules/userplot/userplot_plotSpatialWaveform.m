function S_obj = userplot_plotSpatialWaveform(cav_dir,num,beta,t0)
    S_EcoeffsData = structdata_load(cav_dir,num,'E','coeffs',1,beta); 
    
    [~,t0_ind] = helpers_getClosestMatch(t0,S_EcoeffsData.t);
    S_EcoeffsData.Y = S_EcoeffsData.Y(t0_ind,:);
    S_EcoeffsData.t_ind = t0_ind;
        
    S_obj = plotting_plotData(S_EcoeffsData,'SPCFLD','2Dcoeffs','');
end