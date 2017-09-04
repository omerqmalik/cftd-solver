function userplot_saveMicroFigures(cav_dir,data_dir,num,results_dir,calc_times,pstep,S_coredata)

    %% E-field
    %create directories
    savedir_EfieldFFT  = [results_dir '/EfieldFFT'];
    savedir_EfieldTW   = [results_dir '/EfieldTW'];
    savedir_EcoeffsFFT = [results_dir '/EcoeffsFFT'];
    mkdir(savedir_EfieldFFT);
    mkdir(savedir_EfieldTW);
    mkdir(savedir_EcoeffsFFT);
    
    %Load data
    S_EcoeffsData   = structdata_loadMicro(cav_dir,data_dir,num,'E','coeffs',calc_times,pstep,S_coredata.x0);   %coeffs
    S_EfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'E','field',calc_times,pstep,S_coredata.x0);    %field
    
    %Plot EfieldFFT
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_EfieldFFT,1,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_EfieldFFT,0,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_EfieldFFT,1,1);
    plotting_plotData(S_EcoeffsData,'FFT','2Dfield',savedir_EfieldFFT,0,1);
    
    %plot 2DcoeffsFFT from coeffs
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_EcoeffsFFT,1,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_EcoeffsFFT,1,1);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_EcoeffsFFT,0,0);
    plotting_plotData(S_EcoeffsData,'FFT','2Dcoeffs',savedir_EcoeffsFFT,0,1);
    
    %plot 2DfieldTW from field
    plotting_plotData(S_EfieldData,'TW','2Dfield',savedir_EfieldTW);
    
    
    %% P-field
    %create directories
%     savedir_PfieldFFT  = [results_dir '/PfieldFFT'];
    savedir_PfieldTW   = [results_dir '/PfieldTW'];
%     savedir_PcoeffsFFT = [results_dir '/PcoeffsFFT'];
%     mkdir(savedir_PfieldFFT);
    mkdir(savedir_PfieldTW);
%     mkdir(savedir_PcoeffsFFT);
    
    %Load data
%     S_PcoeffsData   = structdata_loadMicro(cav_dir,data_dir,num,'P','coeffs',calc_times,pstep,S_coredata.x0);   %coeffs
    S_PfieldData    = structdata_loadMicro(cav_dir,data_dir,num,'P','field',calc_times,pstep,S_coredata.x0);    %field
    
%     %Plot EfieldFFT
%     plotting_plotData(S_PcoeffsData,'FFT','2Dfield',savedir_PfieldFFT,1,0);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dfield',savedir_PfieldFFT,0,0);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dfield',savedir_PfieldFFT,1,1);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dfield',savedir_PfieldFFT,0,1);
%     
%     %plot 2DcoeffsFFT from coeffs
%     plotting_plotData(S_PcoeffsData,'FFT','2Dcoeffs',savedir_PcoeffsFFT,1,0);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dcoeffs',savedir_PcoeffsFFT,1,1);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dcoeffs',savedir_PcoeffsFFT,0,0);
%     plotting_plotData(S_PcoeffsData,'FFT','2Dcoeffs',savedir_PcoeffsFFT,0,1);
    
    %plot 2DfieldTW from field
    plotting_plotData(S_PfieldData,'TW','2Dfield',savedir_PfieldTW);
    
    
     %% D-field
%     
%     savedir_DmnAVGABS   = [results_dir '/DmnAVGABS'];
%     savedir_DcoeffsFFT  = [results_dir '/DcoeffsFFT'];
%     
%     mkdir(savedir_DmnAVGABS);
%     mkdir(savedir_DcoeffsFFT);
%     
%     %Load Davgabs data
%     S_DavgabsData   = structdata_loadMicro(cav_dir,data_dir,num,'D','avgabs',calc_times,pstep,S_coredata.x0);
%     
%     %Load Dcoeffs data
%     S_DcoeffsData   = structdata_loadMicro(cav_dir,data_dir,num,'D','coeffs',calc_times,pstep,S_coredata.x0);
%     S_DcoeffsData.Dsave = S_coredata.Dsave;
%     S_DcoeffsData.rframe = 0;
%     
% 
%     %plot DmnAVGABS
%     plotting_plotData(S_DavgabsData,'AVGABS','Dmn',savedir_DmnAVGABS);
%     
%     %plot DcoeffsFFT
%     plotting_plotData(S_DcoeffsData,'FFT','Dcoeffs',savedir_DcoeffsFFT,1,0);
%     plotting_plotData(S_DcoeffsData,'FFT','Dcoeffs',savedir_DcoeffsFFT,0,0);
%     
%     S_DcoeffsArray = diag_getDataArray(S_DcoeffsData);
%     plotting_plotData(S_DcoeffsArray,'FFTmov','Dcoeffs',results_dir,1,0);
% 	plotting_plotData(S_DcoeffsArray,'FFTmov','Dcoeffs',results_dir,0,0);
end