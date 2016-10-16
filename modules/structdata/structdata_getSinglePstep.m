function S_structdata = structdata_getSinglePstep(S_structdata,pstep)
    if strcmp(S_structdata.type,'coeffs')
        S_structdata.Y = S_structdata.Y(:,:,pstep);
    elseif strcmp(S_structdata.type,'field')
        S_structdata.Y = S_structdata.Y(:,pstep);
    end
    S_structdata.calc_times = S_structdata.calc_times(:,pstep);
    S_structdata.pump = S_structdata.pump(pstep);
    
    psteps = 1;
    if length(psteps) > 1
        S_structdata.fname       = [S_structdata.id S_structdata.type '_p' num2str(S_structdata.pump(1)/S_structdata.th) '-' num2str(S_structdata.pump(end)/S_structdata.th) '_t' num2str(S_structdata.t(1)) '-' num2str(S_structdata.t(end))];
    else
        S_structdata.fname       = [S_structdata.id S_structdata.type '_p' num2str(S_structdata.pump(1)/S_structdata.th) '_t' num2str(S_structdata.t(1)) '-' num2str(S_structdata.t(end))];
    end
end