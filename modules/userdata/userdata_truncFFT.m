function [trunc_w,trunc_vec] = userdata_truncFFT(w,fft_w,center,half_width,skip)
    trunc_vec = fft_w(w < (center+half_width) & w > (center-half_width),:);
    trunc_w   = w(w < (center+half_width) & w > (center-half_width));
    
    len = length(trunc_w);
    trunc_vec = trunc_vec(1:skip:len,:);
    trunc_w = trunc_w(1:skip:len);
end