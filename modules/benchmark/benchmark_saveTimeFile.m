function benchmark_saveTimeFile(saveloc,pgroup,calc_times)
    filename = [saveloc '/calctimes_' num2str(pgroup)];
    save(filename,'calc_times');
end