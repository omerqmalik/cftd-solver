function x = helpers_truncVector(x,chunk)
    len_x = size(x,1);
    start = ceil((1-chunk)*len_x)+1;
    stop  = len_x;
    
    x = x(start:stop,:,:);
end