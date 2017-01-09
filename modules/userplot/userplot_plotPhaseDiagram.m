function fldcorr_all = userplot_plotPhaseDiagram(cav_dir,num,teeth,delta)

    %Load coeffs data and also convert to array
    S_EcoeffsAll    = structdata_load(cav_dir,num,'E','coeffs',1);
    [~,CFvecs,~,x] = cavity_loadBasis(S_EcoeffsAll.basis_loc,S_EcoeffsAll.k_a);
    width = cavity_getBasisRange(S_EcoeffsAll.CFvals);
    
    temporal_data   = userdata_calcTemporalField(S_EcoeffsAll.Y,CFvecs,x);
    [w,~,fftw_mag]  = userdata_calcFFT(S_EcoeffsAll.t,temporal_data,S_EcoeffsAll.rframe);
    [w,fftw_mag]    = userdata_truncFFT(w,fftw_mag,S_EcoeffsAll.rframe,width,1);
    
    CFspacing = diff(real(S_EcoeffsAll.CFvals(1:2)));
    
    subwnds = (-CFspacing/2):(CFspacing/delta):(CFspacing/2);
    wnd_ctr = teeth*CFspacing;
    wnd_wdt = repmat(wnd_ctr.',[1,length(subwnds)]) + repmat(subwnds,[length(wnd_ctr),1]);
    
    diff_w = w(2) - w(1);
    fldcorr_all = zeros(size(fftw_mag,2),size(wnd_wdt,1));
    for k = 1:size(fftw_mag,2)
        fld_this = fftw_mag(:,k);
        for i = 1:size(wnd_wdt,1)
            fldcorr  = 0;
            for j = 1:size(wnd_wdt,2)
                wdt_this = wnd_wdt(i,j);
                ind_spacing = fix(wdt_this/diff_w);
                startpts = 1:(length(w)-ind_spacing);
                endpts   = (1+ind_spacing):length(w);
                
                fldstart = fld_this(startpts);
                fldend   = fld_this(endpts);
                fldcorr  = fldcorr + sum(fldstart.*fldend);
            end
            fldcorr_all(k,i) = fldcorr;
        end
    end
end