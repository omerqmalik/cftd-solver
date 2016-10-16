function benchmark_createTimeFiles(saveloc,pgroups,pg_len,pump_sz,len_tvec)
    calc_times = zeros(len_tvec-1,pg_len);
    for i = 1:pgroups
        if i == pgroups
            a = mod(pump_sz,pg_len);
            if a > 0
                calc_times = calc_times(:,1:a);
            end
        end
        if ~(exist(saveloc,'dir') == 7)
            mkdir(saveloc);
        end
        save([saveloc '/calctimes_' num2str(i) '.mat'],'calc_times');
    end
end