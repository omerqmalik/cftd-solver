function pstep_times = benchmark_loadTimeForPstep(saveloc,pstep)
    filename = [saveloc '/psteptimes_' num2str(pstep) '.mat'];
    if exist(filename,'file') == 2
        load(filename,'pstep_times');
    else
        pstep_times = [];
    end
end