function saltplot_laserOutput(cav_dir,num,prel_f,sz)
    [k,a,a2,pump] = salt_getSALTData(cav_dir,num,prel_f,sz);
    pump = pump/min(pump);
    
    figure;
    plot(pump,a2);
    
    figure;
    for i = 1:size(k,2)
        plot(pump(k(:,i)>0),k(k(:,i)>0,i)); hold on;
    end
end