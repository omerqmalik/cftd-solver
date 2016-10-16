function [y_sg] = filter_selectively(x,y,peaks,valleys,isplot)
    if isplot
        plot_peaks_valleys(x,y,peaks,valleys,1);
    end
    
    len_v = size(valleys,1);
    for i = 1:len_v
        if length(y(valleys(i,1):valleys(i,2))) >= 41
            y(valleys(i,1):valleys(i,2)) = sgolayfilt(medfilt1(y(valleys(i,1):valleys(i,2)),10),2,41);
        end
    end
    y_sg = y;
    
    if isplot
        plot_peaks_valleys(x,y,peaks,valleys,0);
    end
end