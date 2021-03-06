function plotting_saveObj(S_obj,results_dir)
    if strcmp(S_obj.type,'2DfieldFFTmov') || strcmp(S_obj.type,'2DcoeffsFFTmov') || strcmp(S_obj.type,'2DfieldTWmov') || strcmp(S_obj.type,'DmnAVGABSmov') || strcmp(S_obj.type,'DcoeffsFFTmov')
        M = S_obj.M;
        save([results_dir '/' S_obj.fname '.mat'],'M');
    elseif strcmp(S_obj.type,'3DfieldFFT') || strcmp(S_obj.type,'2DfieldFFT') || strcmp(S_obj.type,'2DfieldTW') || strcmp(S_obj.type,'2DcoeffsFFT') || strcmp(S_obj.type,'DmnAVGABS') || strcmp(S_obj.type,'DcoeffsFFT') || strcmp(S_obj.type,'3DcoeffsFFT')
        savefig(S_obj.f,[results_dir '/' S_obj.fname '.fig']);
    end
end