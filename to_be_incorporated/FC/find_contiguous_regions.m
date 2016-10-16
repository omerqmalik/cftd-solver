function [peaks,valleys] = find_contiguous_regions(y,thresh_wgt)    %accepts real valued 'y' only
    y = abs(y);
    cutoff = thresh_wgt*median(findpeaks(y));
        
    all_regions = y < cutoff;
    ons_and_offs = diff(all_regions);
    
    val_ons  = sort(find(ons_and_offs == 1))+1;
    val_offs = sort(find(ons_and_offs == -1));
    
    if val_offs(1) < val_ons(1)
        val_ons = [1; val_ons];
    end
    if val_ons(end) > val_offs(end)
        val_offs = [val_offs; length(all_regions)];
    end
    
    g = val_offs(2:end-1) - val_ons(2:end-1) < 100;
    val_ons(logical([0; g; 0])) = [];
    val_offs(logical([0; g; 0])) = [];
    
    peak_ons  = val_offs + 1;
    peak_offs = val_ons - 1;
    
    peak_offs = peak_offs(peak_offs > 0);
    peak_ons  = peak_ons(peak_ons <= length(all_regions));
    
    if peak_offs(1) < peak_ons(1)
        peak_ons = [1; peak_ons];
    end
    if peak_ons(end) > peak_offs(end)
        peak_offs = [peak_offs; length(all_regions)];
    end
    
    peaks = [peak_ons peak_offs];
    valleys = [val_ons val_offs];
    
    peaks(peaks(:,2) - peaks(:,1) < 2,:) = [];
    valleys(valleys(:,2) - valleys(:,1) < 2,:) = [];
end