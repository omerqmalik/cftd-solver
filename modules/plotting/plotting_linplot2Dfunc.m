function [f,h] = plotting_linplot2Dfunc(x,y,S_fig)
    f = plotting_makeFig;
    h = plot(x,y);
    
    xlabel(S_fig.x_label,'interpreter','latex');
    ylabel(S_fig.y_label,'interpreter','latex');
    title(S_fig.title_str,'interpreter','none');
end