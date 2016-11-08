function S_dataArray = structdata_getDataArray(S_structdata)
    len_pump = length(S_structdata.pump);
    S_dataArray(len_pump) = structdata_getPsteps(S_structdata,len_pump);

    for i = 1:(len_pump-1)
        S_dataArray(i) = structdata_getPsteps(S_structdata,i);
    end
end