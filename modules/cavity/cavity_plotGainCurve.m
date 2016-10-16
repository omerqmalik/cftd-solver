function cavity_plotGainCurve(k_a,g_per)
    hold on;
    
    x_lim = get(gca,'xlim');
    w = linspace(x_lim(1),x_lim(2),1000);
    gc = cavity_calcGainCurve(w,k_a,g_per);
    
    box off;
    a_right = axes;
    plot(a_right,w,gc,'color',[.85 .325 .0980]);
    set(gca,'xlim',x_lim);
    set(gca,'ylim',[0 1]);
    
    set(a_right,'YAxisLocation','Right');
    set(a_right,'color','none');
    set(a_right,'xtick',[]);
end
