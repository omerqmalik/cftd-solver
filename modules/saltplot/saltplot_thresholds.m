function saltplot_thresholds(cav_dir,num,prel_f,sz)
    S_coredata = salt_loadParameters(cav_dir,num,prel_f,sz);
    load([S_coredata.data_dir '/rough_scan.mat']);
    load([S_coredata.data_dir '/fine_scan.mat']);
    
    plot(S_rough.k,S_rough.th,'bo'); hold on;
    plot(S_fine.k,S_fine.th,'rx');
    
    yL = get(gca,'ylim');
    for i = 1:length(S_coredata.CFvals)
        line(real(S_coredata.CFvals(i))*[1 1],yL,'color','k','linestyle',':');
    end
end