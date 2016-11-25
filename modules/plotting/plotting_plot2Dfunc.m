function [f,h] = plotting_plot2Dfunc(x,y,S_fig)
    if S_fig.islin == 1
        plot_func = @plot;
    elseif S_fig.islin == 0
        plot_func = @semilogy;
    end

    f = plotting_makeFig;
    
    c = autumn(size(y,2));
    if size(y,2) > 1
        h = plot_func(x,y(:,1),'color',c(1,:));
        hold on;
        for i = 2:size(y,2)
            h = plot_func(x,y(:,i),'color',c(i,:));
        end
    else
        h = plot_func(x,y);
    end
    
    xlabel(S_fig.x_label,'interpreter','latex');
    title(S_fig.title_str,'interpreter','none');
    
    if isfield(S_fig,'y_label')
        ylabel(S_fig.y_label,'interpreter','latex');
    elseif isfield(S_fig,'lgnd')
        legend(S_fig.lgnd,'interpreter','latex');
    end
end