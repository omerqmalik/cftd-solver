function [pk_maxima,pk_locations] = get_FC_peaks(y,peak_ranges)
    len = size(peak_ranges,1);
    pk_maxima = zeros(1,len);
    pk_locations = zeros(1,len);
    
    for i = 1:len
        y_peak = y(peak_ranges(i,1):peak_ranges(i,2));
        
        [all_max,all_loc] = findpeaks(abs(y_peak));
        [~,s] = max(all_max);
        if ~isempty(s)
            pk_maxima(i) = y_peak(all_loc(s));
            pk_locations(i) = peak_ranges(i,1) + all_loc(s) - 1;
        end
    end
    pk_maxima(pk_maxima == 0) = [];
    pk_locations(pk_locations == 0) = [];
end