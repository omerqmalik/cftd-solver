function [w,fft_w,fft_w_mag,fft_max_val,fft_max_loc] = userdata_calcFFT(t,t_sig,rframe)
    L = length(t);
    NFFT = 2^(nextpow2(L)+1);
    
    fft_sig_H = t_sig.*repmat(hanning(length(t_sig)),1,size(t_sig,2),size(t_sig,3));
    fft_w     = flipud(fftshift(fft(fft_sig_H,NFFT,1),1))/L;
    fft_w_mag = 2*abs(fft_w);           %Not really PSD. Factor of 2 due to hanning
    
    t_sample = diff(t(1:2));
    Fs = 1/t_sample;
    v = Fs/2*(-1:2/NFFT:(1-2/NFFT));
    w = 2*pi*v + rframe;
    
    [fft_max_val,loc] = max(fft_w_mag);
    fft_max_loc = w(loc);
end
