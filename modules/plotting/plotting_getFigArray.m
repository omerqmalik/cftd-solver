function S_figArray = plotting_getFigArray(S_dataArray,fig_type,opt_islin,opt_isdec)

    %initialize
    len = length(S_dataArray);
    if strcmp(fig_type,'2DfieldFFT') || strcmp(fig_type,'2DcoeffsFFT') || strcmp(fig_type,'DcoeffsFFT')
        S_figArray(len) = plotting_getFigStruct(S_dataArray(len),fig_type,opt_islin,opt_isdec);
    elseif strcmp(fig_type,'2DfieldTW') || strcmp(fig_type,'DmnAVGABS') || strcmp(fig_type,'2DcoeffsSPCFLD')
        S_figArray(len) = plotting_getFigStruct(S_dataArray(len),fig_type);
    end

    %populate
    for i = 1:(len-1)
        if strcmp(fig_type,'2DfieldFFT') || strcmp(fig_type,'2DcoeffsFFT') || strcmp(fig_type,'DcoeffsFFT')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),fig_type,opt_islin,opt_isdec);
        elseif strcmp(fig_type,'2DfieldTW') || strcmp(fig_type,'DmnAVGABS') || strcmp(fig_type,'2DcoeffsSPCFLD')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),fig_type);
        end
    end
end