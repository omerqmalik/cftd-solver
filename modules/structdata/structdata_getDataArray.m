function S_dataArray = structdata_getDataArray(S_structdata)
    len_pump = length(S_structdata.pump);
    S_dataArray(len_pump) = structdata_getSinglePstep(S_structdata,len_pump);
    S_dataArray(len_pump).fname = S_structdata.fname;

    for i = 1:(len_pump-1)
        S_dataArray(i) = structdata_getSinglePstep(S_structdata,i);
        S_dataArray(i).fname = S_structdata.fname;
    end
end