function t = rawdata_interpT(cT)
    max_T = 0;
    min_Ts = inf;
    
    for i = 1:length(cT)
        if min(cT{i}) > max_T
            max_T = min(cT{i});
        end
        if mean(diff(cT{i})) < min_Ts
            min_Ts = mean(diff(cT{i}));
        end
    end
    
    t = rawdata_getInterpTvec(max_T,max(cT{1}),min_Ts);
end