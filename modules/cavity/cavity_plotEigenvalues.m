function cavity_plotEigenvalues(CFvals,overlay)
    CFvals = real(CFvals);

    if overlay
        x_lim = get(gca,'xlim');
        CFvals = CFvals(CFvals < x_lim(2) & CFvals > x_lim(1));
    end
    
    y_lim = get(gca,'ylim');
    
    for i = 1:length(CFvals)
        line([CFvals(i) CFvals(i)],y_lim,'color','k','linestyle',':');
    end
    set(gca,'ylim',y_lim);
end