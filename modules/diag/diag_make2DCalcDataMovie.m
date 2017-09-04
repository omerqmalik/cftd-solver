function F = diag_make2DCalcDataMovie(dir_name,file_pre)
    files = dir([dir_name '/' file_pre '*.fig']);
    files = rawdata_orderFilesNumerically(files);
    f = openfig([dir_name '/' files(1).name]);
    colormap hot
    
    G = get(gcf,'Position');
    set(gcf,'Position',[G(1) G(2) 2*G(3) 2*G(4)]);
    
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    
    num_frames = length(files);
    F(num_frames) = struct('cdata',[],'colormap',[]);
    drawnow
    F(1) = getframe(gcf);
    close(f);
    
    for i = 2:length(files)
        fn = files(i).name;
        fn = fn(1:(end-4));
        
        f = openfig([dir_name '/' fn '.fig']);
        colormap hot
        G = get(gcf,'Position');
        set(gcf,'Position',[G(1) G(2) 2*G(3) 2*G(4)]);
        
        drawnow
        axis tight manual
        F(i) = getframe(gcf);
        close(f);
    end
    save([dir_name '/Davgabs_mov.mat'],'F');
end