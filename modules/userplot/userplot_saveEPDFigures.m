function userplot_saveEPDFigures(cav_dir,data_dir,num,results_dir,calc_times,inv_time,pstep,S_coredata)

    %create directories
    savedir_2DfieldFFT  = [results_dir '/2DfieldFFT'];
    savedir_2DfieldTW   = [results_dir '/2DfieldTW'];
    
    mkdir(savedir_2DfieldFFT);
    mkdir(savedir_2DfieldTW);
        
    %Load field data
    S_EfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'E','field',calc_times,pstep,S_coredata.x0);
    S_PfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'P','field',calc_times,pstep,S_coredata.x0);
    S_DfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'D','field',inv_time,pstep,S_coredata.x0);
    S_DfieldData.rframe = 0;
    
    %Plot 2DfieldFFT (E)
    plotting_plotData(S_EfieldData,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
    plotting_plotData(S_EfieldData,'FFT','2Dfield',savedir_2DfieldFFT,0,0);
    plotting_plotData(S_EfieldData,'FFT','2Dfield',savedir_2DfieldFFT,1,1);
    plotting_plotData(S_EfieldData,'FFT','2Dfield',savedir_2DfieldFFT,0,1);
    
    %Plot 2DfieldFFT (P)
    plotting_plotData(S_PfieldData,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
    plotting_plotData(S_PfieldData,'FFT','2Dfield',savedir_2DfieldFFT,0,0);
    plotting_plotData(S_PfieldData,'FFT','2Dfield',savedir_2DfieldFFT,1,1);
    plotting_plotData(S_PfieldData,'FFT','2Dfield',savedir_2DfieldFFT,0,1);
    
    %Plot 2DfieldFFT (D)
    plotting_plotData(S_DfieldData,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
    plotting_plotData(S_DfieldData,'FFT','2Dfield',savedir_2DfieldFFT,0,0);

    %plot 2DfieldTW (E,P,D)
    plotting_plotData(S_EfieldData,'TW','2Dfield',savedir_2DfieldTW);
    plotting_plotData(S_PfieldData,'TW','2Dfield',savedir_2DfieldTW);
    plotting_plotData(S_DfieldData,'TW','2Dfield',savedir_2DfieldTW);
end