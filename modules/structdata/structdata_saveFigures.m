function structdata_saveFigures(cav_dir,num,results_dir)

    %plot 3Dfft
    S_structdata  = structdata_load(cav_dir,num,'E','coeffs',1);
    S_fig   = plotting_getFigStruct(S_structdata,'fft_field');
    S_fig.f = structdata_plotFFT(S_structdata,S_fig);
    plotting_saveObj(S_fig,results_dir);
    close(S_fig.f);
    
    %plot 2Dffts (lin and log)
    pump = S_structdata.pump;
    saveloc_2D = [results_dir '/2Dfft'];
    saveloc_tw = [results_dir '/tw'];
    if ~(exist(saveloc_2D,'dir') == 7)
        mkdir(saveloc_2D);
    end
    if ~(exist(saveloc_tw,'dir') == 7)
        mkdir(saveloc_tw);
    end
    for i = 1:length(pump)
        S_thisdata = structdata_getSinglePstep(S_structdata,i);

        S_fig      = plotting_getFigStruct(S_thisdata,'fft_field',1,0);
        S_fig.f    = structdata_plotFFT(S_thisdata,S_fig);
        plotting_saveObj(S_fig,saveloc_2D);
        close(S_fig.f);
        
        S_fig      = plotting_getFigStruct(S_thisdata,'fft_field',0,0);
        S_fig.f    = structdata_plotFFT(S_thisdata,S_fig);
        plotting_saveObj(S_fig,saveloc_2D);
        close(S_fig.f);
        
        S_fig      = plotting_getFigStruct(S_thisdata,'fft_field',1,1);
        S_fig.f    = structdata_plotFFT(S_thisdata,S_fig);
        plotting_saveObj(S_fig,saveloc_2D);
        close(S_fig.f);
        
        S_fig      = plotting_getFigStruct(S_thisdata,'fft_field',0,1);
        S_fig.f    = structdata_plotFFT(S_thisdata,S_fig);
        plotting_saveObj(S_fig,saveloc_2D);
        close(S_fig.f);
        
        S_fig      = plotting_getFigStruct(S_thisdata,'tw');
        S_fig.f    = structdata_plotTemporalWaveform(S_thisdata,S_fig);
        plotting_saveObj(S_fig,saveloc_tw);
        close(S_fig.f);
    end
    
    %plot 2Dfft movie (lin and log)
    S_dataArray = structdata_getDataArray(S_structdata);
    S_figArray  = plotting_getFigArray(S_dataArray,'fft',1,1);
    S_movie     = structdata_makeMovie('fft',S_dataArray,S_figArray);
    plotting_saveObj(S_movie,results_dir);
    
    S_figArray  = plotting_getFigArray(S_dataArray,'fft',0,1);
    S_movie     = structdata_makeMovie('fft',S_dataArray,S_figArray);
    plotting_saveObj(S_movie,results_dir);
    
%     %plot E_t
%     S_fig   = plotting_getFigStruct(S_structdata,'tw');
%     S_fig.f = structdata_plotTemporalWaveform(S_structdata,S_fig);
%     plotting_saveObj(S_fig,results_dir);
%     close(S_fig.f);
    
    %plot E_t movie
    S_figArray  = plotting_getFigArray(S_dataArray,'tw');
    S_movie     = structdata_makeMovie('tw',S_dataArray,S_figArray);
    plotting_saveObj(S_movie,results_dir);
end