function CF_range = cavity_getBasisRange(CFvals)
    CFvals = real(CFvals);

    width = abs(max(CFvals)-min(CFvals))/2;
    CF_range =  width + 0.1*width;  %Add a fudge factor
end
