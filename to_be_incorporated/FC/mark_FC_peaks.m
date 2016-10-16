function mark_FC_peaks(x,FC_pkmax,FC_pkloc)
    hold on;
    for i = 1:length(FC_pkloc)
        plot(x(FC_pkloc(i)),FC_pkmax(i),'color','r','marker','x');
    end
end