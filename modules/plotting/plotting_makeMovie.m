function F = plotting_makeMovie(plot_func,frame_var1,frame_var2)
    f = plot_func(frame_var1(1),frame_var2(1));
    
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    
    num_frames = length(frame_var1);
    F(num_frames) = struct('cdata',[],'colormap',[]);
    drawnow
    
    %uncomment line below to control xlim/ylim
%     set(gca,'ylim',[0 50]);
    
    F(1) = getframe(gcf);
    close(f);
    for j = 2:num_frames
        f = plot_func(frame_var1(j),frame_var2(j));
        axis tight manual
        drawnow
        
        %uncomment line below to control xlim/ylim
%         set(gca,'ylim',[0 50]);

        F(j) = getframe(gcf);
        close(f);
    end
end