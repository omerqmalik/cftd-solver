function calc_times = benchmark_loadTimeFile(saveloc,pgroup)
    load([saveloc,'/calctimes_',num2str(pgroup),'.mat'],'calc_times');
end