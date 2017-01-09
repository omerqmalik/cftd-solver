function userplot_plot3DFFT(cav_dir,num)

    %Load coeffs data and also convert to array
    S_EcoeffsAll   = structdata_load(cav_dir,num,'E','coeffs',1);
        
    %Plot 3DfieldFFT
    plotting_plotData(S_EcoeffsAll,'FFT','3Dfield','');
end