function [x,y] = rawdata_removeDuplicates(x,y)
    [x,s] = unique(x);
    y = y(s);
end