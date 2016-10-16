function benchmark_consolidateTimes(times_dir)
    timefiles = dir([times_dir '/calctimes_*.mat']);
    timefiles = rawdata_orderFilesNumerically(timefiles);

    times = [];
    for i = 1:length(timefiles)
        load([times_dir '/' timefiles(i).name],'calc_times');
        times = [times calc_times];
    end
    calc_times = times;
    save([times_dir '/calctimes.mat'],'calc_times');
end