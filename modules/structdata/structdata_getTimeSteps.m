function S_structdata = structdata_getTimeSteps(S_structdata,time_ind)
    S_structdata.Y = S_structdata.Y(time_ind,:);
    S_structdata.t_ind = time_ind;
end