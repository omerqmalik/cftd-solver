function [f,h] = plotting_plotMatrix(Y,S_fig)
    f = plotting_makeFig;
    
    h = imagesc(Y);
    colormap hot;
    
    xlabel(S_fig.x_label,'interpreter','latex');
    ylabel(S_fig.y_label,'interpreter','latex');
    title(S_fig.title_str,'interpreter','none');
    
    h1 = colorbar;
    ylabel(h1,S_fig.z_label,'interpreter','latex');
end