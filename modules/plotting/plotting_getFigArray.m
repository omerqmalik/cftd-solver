function S_figArray = plotting_getFigArray(S_dataArray,fig_type,opt_islin,opt_isdec)

    %initialize
    len_pump = length(S_dataArray);
    if strcmp(fig_type,'2DfieldFFT') || strcmp(fig_type,'2DcoeffsFFT')
        S_figArray(len_pump) = plotting_getFigStruct(S_dataArray(len_pump),fig_type,opt_islin,opt_isdec);
    elseif strcmp(fig_type,'2DfieldTW')
        S_figArray(len_pump) = plotting_getFigStruct(S_dataArray(len_pump),fig_type);
    end

    %populate
    for i = 1:(len_pump-1)
        if strcmp(fig_type,'2DfieldFFT') || strcmp(fig_type,'2DcoeffsFFT')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),fig_type,opt_islin,opt_isdec);
        elseif strcmp(fig_type,'2DfieldTW')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),fig_type);
        end
    end
end