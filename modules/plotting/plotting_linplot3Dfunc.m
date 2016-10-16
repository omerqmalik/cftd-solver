function [f,h] = plotting_linplot3Dfunc(x,y,Z,S_fig)
    f = plotting_makeFig;
    h = surf(x,y,Z);
    shading interp;
    view(90,-90);
    
    xlabel(S_fig.x_label,'interpreter','latex');
    ylabel(S_fig.y_label,'interpreter','latex');
    
    title(S_fig.title_str,'interpreter','none');
    
    h = colorbar;
    ylabel(h,S_fig.z_label);
end