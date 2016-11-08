function structdata_saveFigures(cav_dir,num,results_dir)

    %create directories
    savedir_2DfieldFFT  = [results_dir '/2DfieldFFT'];
    savedir_2DfieldTW   = [results_dir '/2DfieldTW'];
    savedir_2DcoeffsFFT = [results_dir '/2DcoeffsFFT'];
    if ~(exist(savedir_2DfieldFFT,'dir') == 7)
        mkdir(savedir_2DfieldFFT);
    end
    if ~(exist(savedir_2DfieldTW,'dir') == 7)
        mkdir(savedir_2DfieldTW);
    end
    if ~(exist(savedir_2DcoeffsFFT,'dir') == 7)
        mkdir(savedir_2DcoeffsFFT);
    end

    %Load coeffs data and also convert to array
    S_structdataAll = structdata_load(cav_dir,num,'E','coeffs',1);
    S_dataArray     = structdata_getDataArray(S_structdataAll);
    
    %Plot 3DfieldFFT
    plotting_plotData(S_structdataAll,'FFT','3Dfield',results_dir);
    
    pump = S_structdataAll.pump;
    for i = 1:length(pump)
        %Get individual pump step structs
        S_structdataSingle = structdata_getPsteps(S_structdataAll,i);
 
        %plot 2DfieldFFT (lin and log)
        plotting_plotData(S_structdataSingle,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
        plotting_plotData(S_structdataSingle,'FFT','2Dfield',savedir_2DfieldFFT,0,0);
        plotting_plotData(S_structdataSingle,'FFT','2Dfield',savedir_2DfieldFFT,1,1);
        plotting_plotData(S_structdataSingle,'FFT','2Dfield',savedir_2DfieldFFT,0,1);
        
        %plot 2DcoeffsFFT (log)
        plotting_plotData(S_structdataSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,0);
        plotting_plotData(S_structdataSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,1);

        %plot 2DfieldTW (from coeffs)
        plotting_plotData(S_structdataSingle,'TW','2Dfield',savedir_2DfieldTW);
    end
    
    %plot 2DfieldFFTmov (lin and log)
    plotting_plotData(S_dataArray,'FFTmov','2Dfield',results_dir,1,1)
    plotting_plotData(S_dataArray,'FFTmov','2Dfield',results_dir,0,1)
    
    %plot 2DcoeffsFFTmov (log)
    plotting_plotData(S_dataArray,'FFTmov','2Dcoeffs',results_dir,0,1);
        
    %plot 2DfieldTWmov (from coeffs)
    plotting_plotData(S_dataArray,'TWmov','2Dfield',results_dir);
    
    %plot 2DfieldTWmov (from field)
    S_structdataAll = structdata_load(cav_dir,num,'E','field',1);
    S_dataArray  = structdata_getDataArray(S_structdataAll);
    plotting_plotData(S_dataArray,'TWmov','2Dfield',results_dir);
    
    %plot 2DfieldTW (from field)
    pump = S_structdataAll.pump;
    for i = 1:length(pump)
        S_structdataSingle = structdata_getPsteps(S_structdataAll,i);
        
        plotting_plotData(S_structdataSingle,'TW','2Dfield',savedir_2DfieldTW);
    end
end