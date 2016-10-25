function benchmark_plotPstepProgress()
    %load pump and calcs details
    S_setupdata = setup_loadParameters('.');
    len = length(S_setupdata);
    
    %load psteps being calculated
    timefiles = rawdata_orderFilesNumerically(dir([S_setupdata(1).times_dir '/psteptimes_*.mat']));
    load([S_setupdata(1).times_dir '/' timefiles(1).name],'pstep_times');
    
    %create empty array to hold all times
    n_rows    = length(S_setupdata(1).pump)*length(pstep_times);
    n_columns = len;
    chunk_len = length(pstep_times);
    times_mat = zeros(n_rows,n_columns);
    
    for i = 1:len
        [timefiles,fnum] = rawdata_orderFilesNumerically(dir([S_setupdata(i).times_dir '/psteptimes_*.mat']));
        for j = 1:length(timefiles)
            load([S_setupdata(i).times_dir '/' timefiles(j).name],'pstep_times');
            times_mat((fnum(j)-1)*chunk_len+1:chunk_len*fnum(j),i) = pstep_times;
        end
    end
    total_time = sum(sum(times_mat))/60/60;
    
    times_mat(:,end+1) = zeros(size(times_mat,1),1);
    times_mat(end+1,:) = zeros(1,size(times_mat,2));
    lenx = size(times_mat,2);
    leny = size(times_mat,1);

    surf(1:lenx,1:leny,times_mat);
    shading flat;
        
    view(0,90);
    set(gca,'xtick',1:lenx);
    set(gca,'ytick',[1 chunk_len*(1:length(S_setupdata(1).pump))]);
    set(gca,'yticklabel',[0 1:length(S_setupdata(1).pump)]);
    xlabel('calc #');
    ylabel('pstep');
    
    h = colorbar;
    ylabel(h,'time (s)');
    
    p = pwd;
    j = find(p == '/');
    title([ 'Calculation Progress (' p(j(length(j)-1)+1:j(length(j))-1) ' - ' p(j(length(j))+1:end) ') - ' num2str(total_time) 'h'],'interpreter','none');
end