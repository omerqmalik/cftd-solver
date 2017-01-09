function S_obj = plotting_plotData(S_structdata,func_type,func_name,save_dir,opt_islin,opt_isdec)
    if strcmp(func_type,'FFT')
        if strcmp(func_name,'3Dfield') || strcmp(func_name,'3Dcoeffs')
            S_obj = plotting_getFigStruct(S_structdata,[func_name 'FFT']);
        elseif strcmp(func_name,'2Dfield') || strcmp(func_name,'2Dcoeffs') || strcmp(func_name,'Dcoeffs')
            S_obj = plotting_getFigStruct(S_structdata,[func_name 'FFT'],opt_islin,opt_isdec);
        end
        S_obj.f = structdata_plotFFT(S_structdata,S_obj);
    elseif strcmp(func_type,'TW')
        if strcmp(func_name,'2Dfield')
            S_obj = plotting_getFigStruct(S_structdata,[func_name 'TW']);
        end
        S_obj.f = structdata_plotTemporalWaveform(S_structdata,S_obj);
    elseif strcmp(func_type,'SPCFLD')
        if strcmp(func_name,'2Dcoeffs')
            S_obj = plotting_getFigStruct(S_structdata,[func_name 'SPCFLD']);
        end
        S_obj.f = structdata_plotSpatialWaveform(S_structdata,S_obj);
    elseif strcmp(func_type,'FFTmov')
        S_figArray = plotting_getFigArray(S_structdata,[func_name 'FFT'],opt_islin,opt_isdec);
        S_obj      = plotting_getMovStruct([func_name 'FFTmov'],S_structdata,opt_islin,opt_isdec);
        S_obj.M    = plotting_makeMovie(@(S_thisdata,S_thisfig) structdata_plotFFT(S_thisdata,S_thisfig),S_structdata,S_figArray);
    elseif strcmp(func_type,'TWmov')
        S_figArray = plotting_getFigArray(S_structdata,[func_name 'TW']);
        S_obj      = plotting_getMovStruct([func_name 'TWmov'],S_structdata);
        S_obj.M    = plotting_makeMovie(@(S_thisdata,S_thisfig) structdata_plotTemporalWaveform(S_thisdata,S_thisfig),S_structdata,S_figArray);
    elseif strcmp(func_type,'AVGABS')
        if strcmp(func_name,'Dmn')
            S_obj = plotting_getFigStruct(S_structdata,[func_name 'AVGABS']);
        end
        S_obj.f = diag_plot2DAvgAbs(S_structdata,S_obj);
    elseif strcmp(func_type,'AVGABSmov')
        S_figArray = plotting_getFigArray(S_structdata,[func_name 'AVGABS']);
        S_obj      = plotting_getMovStruct([func_name 'AVGABSmov'],S_structdata);
        S_obj.M    = plotting_makeMovie(@(S_thisdata,S_thisfig) diag_plot2DAvgAbs(S_thisdata,S_thisfig),S_structdata,S_figArray);
    elseif strcmp(func_type,'SPCFLDmov')
        S_figArray  = plotting_getFigArray(S_structdata,'2DcoeffsSPCFLD');
        S_obj       = plotting_getMovStruct('2DcoeffsSPCFLDmov',S_structdata);
        S_obj.M     = plotting_makeMovie(@(S_thisdata,S_thisfig) structdata_plotSpatialWaveform(S_thisdata,S_thisfig),S_structdata,S_figArray);
    end
    
    if ~strcmp(save_dir,'')
        plotting_saveObj(S_obj,save_dir);
        
        if strcmp(func_type,'FFT') || strcmp(func_type,'TW') || strcmp(func_type,'AVGABS')
            close(S_obj.f);
        end
    end
end