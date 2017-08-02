function time_inds = helpers_getSptlFldMovTimeIndices(S_EcoeffsAll,t0,t1,numframes)
    t = S_EcoeffsAll.t;
    
    [~,t0_ind] = helpers_getClosestMatch(t0,t);
    [~,t1_ind] = helpers_getClosestMatch(t1,t);
    
    d = fix((t1_ind - t0_ind)/(numframes - 1));
    
    if d > 0
        time_inds = t0_ind:d:t1_ind;
        if time_inds(end) < t1_ind
            time_inds = [time_inds (time_inds(end)+d)];
        end
    else
        time_inds = t0_ind;
    end
end