function plot_before_after(x,y,y_sg,peaks,valleys,FC_pkloc,FC_pkmax)
    figure;
    for i = 1:size(peaks,1)
        x_this = x(peaks(i,1):peaks(i,2));
        y_this = y(peaks(i,1):peaks(i,2));
        ysg_this = y_sg(peaks(i,1):peaks(i,2));
        plot(x_this,ysg_this,'b-'); hold on;
        plot(x_this,y_this,'b.');
    end
    
    for i = 1:size(valleys,1)
        x_this = x(valleys(i,1):valleys(i,2));
        y_this = y(valleys(i,1):valleys(i,2));
        ysg_this = y_sg(valleys(i,1):valleys(i,2));
        plot(x_this,ysg_this,'g-'); hold on;
        plot(x_this,y_this,'r.');
    end
    
    mark_FC_peaks(x,FC_pkloc,FC_pkmax);
end