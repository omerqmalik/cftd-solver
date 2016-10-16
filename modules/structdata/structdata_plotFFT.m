function f = structdata_plotFFT(S_structdata,S_fig)
    if strcmp(S_structdata.type,'coeffs')
        [~,CFvecs] = cavity_loadBasis(S_structdata.basis_loc,S_structdata.k_a);
        field_t    = userdata_calcTemporalField(S_structdata.Y,CFvecs);
    elseif strcmp(S_structdata.type,'field')
        field_t = S_structdata.Y;
    end
    [w,~,fftw_mag] = userdata_calcFFT(S_structdata.t,field_t,S_structdata.rframe);
    [w,fftw_mag]   = userdata_truncFFT(w,fftw_mag,S_structdata.rframe,cavity_getBasisRange(S_structdata.CFvals),1);
    
    if size(fftw_mag,2) > 1
        f = plotting_linplot3Dfunc(S_structdata.pump/S_structdata.th,w,fftw_mag,S_fig);
    else
        if S_fig.islin
            f = plotting_linplot2Dfunc(w,fftw_mag,S_fig);
        else
            f = plotting_logplot2Dfunc(w,fftw_mag,S_fig);
        end
        if S_fig.isdec
            cavity_plotEigenvalues(S_structdata.CFvals,1);
            cavity_plotGainCurve(S_structdata.k_a,S_structdata.g_per);
        end
    end
end