function S_dataArray = diag_getDataArray(S_structdata)
    S_dataArray(length(S_structdata.Dsave)) = S_structdata;
    for i = 1:length(S_structdata.Dsave)
        S_dataArray(i) = S_structdata;
        S_dataArray(i).Y = S_structdata.Y(:,i);
        S_dataArray(i).Dsave = S_structdata.Dsave(i);
    end
end