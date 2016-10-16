function S_figArray = plotting_getFigArray(S_dataArray,fig_type,varargin)
    len_pump = length(S_dataArray);
    if strcmp(fig_type,'fft')
        S_figArray(len_pump) = plotting_getFigStruct(S_dataArray(len_pump),'fft_field',varargin{1},varargin{2});
    elseif strcmp(fig_type,'tw')
        S_figArray(len_pump) = plotting_getFigStruct(S_dataArray(len_pump),'tw');
    end

    for i = 1:(len_pump-1)
        if strcmp(fig_type,'fft')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),'fft_field',varargin{1},varargin{2});
        elseif strcmp(fig_type,'tw')
            S_figArray(i) = plotting_getFigStruct(S_dataArray(i),'tw');
        end
    end
end