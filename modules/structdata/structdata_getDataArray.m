function S_dataArray = structdata_getDataArray(S_structdata,type,varargin)
    if strcmp(type,'pump')
        len_pump = length(S_structdata.pump);
        S_dataArray(len_pump) = structdata_getPumpSteps(S_structdata,len_pump);

        for i = 1:(len_pump-1)
            S_dataArray(i) = structdata_getPumpSteps(S_structdata,i);
        end
    elseif strcmp(type,'time')
        time_inds = varargin{1};
        len_times = length(time_inds);
        S_dataArray(len_times) = structdata_getTimeSteps(S_structdata,time_inds(end));
        
        for i = 1:(len_times - 1)
            S_dataArray(i) = structdata_getTimeSteps(S_structdata,time_inds(i));
        end
    end
end