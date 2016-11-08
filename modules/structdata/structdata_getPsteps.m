function S_structdata = structdata_getPsteps(S_structdata,psteps)
    if strcmp(S_structdata.type,'coeffs')
        S_structdata.Y = S_structdata.Y(:,:,psteps);
    elseif strcmp(S_structdata.type,'field')
        S_structdata.Y = S_structdata.Y(:,psteps);
    end
    S_structdata.calc_times = S_structdata.calc_times(:,psteps);
    S_structdata.pump = S_structdata.pump(psteps);
end