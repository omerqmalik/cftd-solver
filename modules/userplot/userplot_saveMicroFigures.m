function userplot_saveMicroFigures(cav_dir,data_dir,num,results_dir,calc_times,pstep,S_coredata)

    %create directories
    savedir_2DfieldFFT  = [results_dir '/2DfieldFFT'];
    savedir_2DfieldTW   = [results_dir '/2DfieldTW'];
    savedir_2DcoeffsFFT = [results_dir '/2DcoeffsFFT'];
    savedir_DmnAVGABS   = [results_dir '/DmnAVGABS'];
    savedir_DcoeffsFFT  = [results_dir '/DcoeffsFFT'];
    
    mkdir(savedir_2DfieldFFT);
    mkdir(savedir_2DfieldTW);
    mkdir(savedir_2DcoeffsFFT);
    mkdir(savedir_DmnAVGABS);
    mkdir(savedir_DcoeffsFFT);
    
    %Load coeffs data
    S_EcoeffsData   = structdata_loadMicro(cav_dir,data_dir,num,'E','coeffs',calc_times,pstep,S_coredata.x0);
    
    %Load field data
    S_EfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'E','field',calc_times,pstep,S_coredata.x0);
    
    %Load Davgabs data
    S_DavgabsData   = structdata_loadMicro(cav_dir,data_dir,num,'D','avgabs',calc_times,pstep,S_coredata.x0);
    
    %Load Dcoeffs data
    S_DcoeffsData   = structdata_loadMicro(cav_dir,data_dir,num,'D','coeffs',calc_times,pstep,S_coredata.x0);
    S_DcoeffsData.Dsave = S_coredata.Dsave;
    S_DcoeffsData.rframe = 0;
    
    %Plot 2DfieldFFT
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_2DfieldFFT,0,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_2DfieldFFT,1,1);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_2DfieldFFT,0,1);

    %plot 2DcoeffsFFT from coeffs (log)
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,1,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,1,1);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,1);

    %plot 2DfieldTW from field
    plotting_plotData(S_EfieldData,'TW','2Dfield',savedir_2DfieldTW);

    %plot DmnAVGABS
    plotting_plotData(S_DavgabsData,'AVGABS','Dmn',savedir_DmnAVGABS);
    
    %plot DcoeffsFFT
    plotting_plotData(S_DcoeffsData,'FFT','Dcoeffs',savedir_DcoeffsFFT,1,0);
    plotting_plotData(S_DcoeffsData,'FFT','Dcoeffs',savedir_DcoeffsFFT,0,0);
    
    S_DcoeffsArray = diag_getDataArray(S_DcoeffsData);
    plotting_plotData(S_DcoeffsArray,'FFTmov','Dcoeffs',results_dir,1,0);
	plotting_plotData(S_DcoeffsArray,'FFTmov','Dcoeffs',results_dir,0,0);
end