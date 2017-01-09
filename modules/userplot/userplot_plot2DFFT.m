function userplot_plot2DFFT(cav_dir,num,beta,islin,isdec)

    %Load coeffs data and also convert to array
    S_EcoeffsAll   = structdata_load(cav_dir,num,'E','coeffs',1,beta);
    
    %plot 2DfieldFFT from coeffs (lin and log)
    plotting_plotData(S_EcoeffsAll,'FFT','2Dfield','',islin,isdec);
end