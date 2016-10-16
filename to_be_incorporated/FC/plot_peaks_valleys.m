function plot_peaks_valleys(x,y,peaks,valleys,is_new)
    if(is_new)
        figure;
        c = 'b';
    else
        c = 'r';
    end
    subplot(3,1,1);
    plot(x,y,'color',c);
    hold on;
    
    subplot(3,1,2);
    for i = 1:size(peaks,1)
        x_this = x(peaks(i,1):peaks(i,2));
        y_this = y(peaks(i,1):peaks(i,2));
        plot(x_this,y_this,'color',c);
        hold on;
    end
    
    subplot(3,1,3);
    for i = 1:size(valleys,1)
        x_this = x(valleys(i,1):valleys(i,2));
        y_this = y(valleys(i,1):valleys(i,2));
        plot(x_this,y_this,'color',c);
        hold on;
    end
end