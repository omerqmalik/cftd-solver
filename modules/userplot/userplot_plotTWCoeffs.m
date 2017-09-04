function userplot_plotTWCoeffs(cav_dir,num,beta)

    %Load coeffs data and also convert to array
    S_EcoeffsAll   = structdata_load(cav_dir,num,'E','coeffs',1,0,beta);
    
    %plot 2DfieldFFT from coeffs (lin and log)
    plotting_plotData(S_EcoeffsAll,'TW','2Dcoeffs','',1,0);
end