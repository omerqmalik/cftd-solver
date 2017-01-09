function S_structdata = structdata_getPumpSteps(S_structdata,psteps)
    if strcmp(S_structdata.type,'coeffs') || strcmp(S_structdata.type,'avgabs')
        S_structdata.Y = S_structdata.Y(:,:,psteps);
    elseif strcmp(S_structdata.type,'field')
        S_structdata.Y = S_structdata.Y(:,psteps);
    end
    S_structdata.calc_times = S_structdata.calc_times(:,psteps);
    S_structdata.pump = S_structdata.pump(psteps);
    S_structdata.psteps = psteps;
end