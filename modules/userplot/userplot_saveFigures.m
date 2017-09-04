function userplot_saveFigures(cav_dir,num,x0,results_dir)

    %create directories
    savedir_2DfieldFFT  = [results_dir '/2DfieldFFT'];
    savedir_2DfieldTW   = [results_dir '/2DfieldTW'];
    savedir_2DcoeffsFFT = [results_dir '/2DcoeffsFFT'];
    savedir_DmnAVGABS   = [results_dir '/DmnAVGABS'];
    if ~(exist(savedir_2DfieldFFT,'dir') == 7)
        mkdir(savedir_2DfieldFFT);
    end
    if ~(exist(savedir_2DfieldTW,'dir') == 7)
        mkdir(savedir_2DfieldTW);
    end
    if ~(exist(savedir_2DcoeffsFFT,'dir') == 7)
        mkdir(savedir_2DcoeffsFFT);
    end
%     if ~(exist(savedir_DmnAVGABS,'dir') == 7)
%         mkdir(savedir_DmnAVGABS);
%     end

    %Load coeffs data and also convert to array
    S_EcoeffsAll   = structdata_load(cav_dir,num,'E','coeffs',1,x0);
    S_EcoeffsArray = structdata_getDataArray(S_EcoeffsAll,'pump');
    
    %Load field data and also convert to array
    S_EfieldAll    = structdata_load(cav_dir,num,'E','field',1,x0);
    S_EfieldArray  = structdata_getDataArray(S_EfieldAll,'pump');
    
%     %Load Davgabs data and also convert to array
%     S_DavgabsAll   = structdata_load(cav_dir,num,'D','avgabs',0,x0);
%     S_DavgabsArray = structdata_getDataArray(S_DavgabsAll,'pump');
    
    if strcmp(S_EcoeffsAll.pump_type,'hysteresis')
        pump_mid = ceil(length(S_EcoeffsAll.pump)/2);
        pump_up = 1:pump_mid;
        pump_down = pump_mid:length(S_EcoeffsAll.pump);
        S_EcoeffsUp = structdata_getPumpSteps(S_EcoeffsAll,pump_up);
        S_EcoeffsDown = structdata_getPumpSteps(S_EcoeffsAll,pump_down);
        
        %Plot 3DfieldFFT
        plotting_plotData(S_EcoeffsUp,'FFT','3Dfield',results_dir);
        plotting_plotData(S_EcoeffsDown,'FFT','3Dfield',results_dir);

        %Plot 3DcoeffsFFT
        plotting_plotData(S_EcoeffsUp,'FFT','3Dcoeffs',results_dir);
        plotting_plotData(S_EcoeffsDown,'FFT','3Dcoeffs',results_dir);
    elseif strcmp(S_EcoeffsAll.pump_type,'simple')
        plotting_plotData(S_EcoeffsAll,'FFT','3Dfield',results_dir);
        plotting_plotData(S_EcoeffsAll,'FFT','3Dcoeffs',results_dir);
    end
    
    %plot 2DfieldFFTmov from coeffs (lin and log)
    plotting_plotData(S_EcoeffsArray,'FFTmov','2Dfield',results_dir,1,1);
    plotting_plotData(S_EcoeffsArray,'FFTmov','2Dfield',results_dir,0,1);
    
    %plot 2DcoeffsFFTmov from coeffs (log)
    plotting_plotData(S_EcoeffsArray,'FFTmov','2Dcoeffs',results_dir,0,1);
    plotting_plotData(S_EcoeffsArray,'FFTmov','2Dcoeffs',results_dir,1,1);
    
    %plot 2DfieldTWmov from field
    plotting_plotData(S_EfieldArray,'TWmov','2Dfield',results_dir);
        
%     %plot DmnAVGABSmov
%     plotting_plotData(S_DavgabsArray,'AVGABSmov','Dmn',results_dir);
    
    pump = S_EcoeffsAll.pump;
    for i = 1:length(pump)
        %Get individual pump step structs
        S_EcoeffsSingle = structdata_getPumpSteps(S_EcoeffsAll,i);
        S_EfieldSingle  = structdata_getPumpSteps(S_EfieldAll,i);
%         S_DavgabsSingle = structdata_getPumpSteps(S_DavgabsAll,i);
 
        %plot 2DfieldFFT from coeffs (lin and log)
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dfield',savedir_2DfieldFFT,1,0);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dfield',savedir_2DfieldFFT,0,0);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dfield',savedir_2DfieldFFT,1,1);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dfield',savedir_2DfieldFFT,0,1);
        
        %plot 2DcoeffsFFT from coeffs (log)
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,1,0);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,1,1);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,0);
        plotting_plotData(S_EcoeffsSingle,'FFT','2Dcoeffs',savedir_2DcoeffsFFT,0,1);
        
        %plot 2DfieldTW from field
        plotting_plotData(S_EfieldSingle,'TW','2Dfield',savedir_2DfieldTW);
        
%         %plot DmnAVGABS
%         plotting_plotData(S_DavgabsSingle,'AVGABS','Dmn',savedir_DmnAVGABS);
    end
end