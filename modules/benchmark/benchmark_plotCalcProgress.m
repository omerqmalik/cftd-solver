function benchmark_plotCalcProgress()
    S_setupdata = setup_loadParameters('.');
    
    len = length(S_setupdata);
    times_vec = zeros(length(S_setupdata(1).pump),len);
    for i = 1:len
        timefiles = dir([S_setupdata(i).times_dir '/calctimes_*.mat']);
        timefiles = rawdata_orderFilesNumerically(timefiles);
        
        start = 1;
        for ii = 1:length(timefiles)
            load([S_setupdata(i).times_dir '/' timefiles(ii).name],'calc_times');
            stop = start + size(calc_times,2) - 1;
            
            times_vec(start:stop,i) = sum(calc_times,1).';
            start = stop + 1;
        end
    end
    total_time = sum(sum(times_vec))/60/60;
    
    y = S_setupdata(1).pump/S_setupdata(1).th;
    if length(size(times_vec)) == 2 && size(times_vec,1) == 1
        len = size(times_vec,2);
        plot(1:len,times_vec/60);
        xlabel('calc #');
        ylabel('time (minutes)');
    else
        times_vec(:,i+1) = zeros(length(S_setupdata(1).pump),1);
        len = size(times_vec,2);

        surf(1:len,y,times_vec/60);
        
        view(0,90);
        set(gca,'xlim',[1,len]);
        set(gca,'ylim',[min(y), max(y)]);
        xlabel('calc #');
        ylabel({'$\beta$'},'interpreter','latex');
        
        h = colorbar;
        ylabel(h,'time (minutes)');
    end
    
    p = pwd;
    j = find(p == '/');
    title([ 'Calculation Progress (' p(j(length(j)-1)+1:j(length(j))-1) ' - ' p(j(length(j))+1:end) ') - ' num2str(total_time) 'h'],'interpreter','none');
end
