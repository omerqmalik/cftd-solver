function f = structdata_plotFFT(S_structdata,S_fig)
    if strcmp(S_fig.type,'2DfieldFFT') || strcmp(S_fig.type,'3DfieldFFT')
        if strcmp(S_structdata.type,'coeffs')
            [~,CFvecs]    = cavity_loadBasis(S_structdata.basis_loc,S_structdata.k_a);
            temporal_data = userdata_calcTemporalField(S_structdata.Y,CFvecs);
        elseif strcmp(S_structdata.type,'field')
            temporal_data = S_structdata.Y;
        end
    elseif strcmp(S_fig.type,'2DcoeffsFFT')
        temporal_data = S_structdata.Y;
    end
        
    [w,~,fftw_mag] = userdata_calcFFT(S_structdata.t,temporal_data,S_structdata.rframe);
    [w,fftw_mag]   = userdata_truncFFT(w,fftw_mag,S_structdata.rframe,cavity_getBasisRange(S_structdata.CFvals),1);
    
    if strcmp(S_fig.type,'3DfieldFFT')
        f = plotting_linplot3Dfunc(S_structdata.pump/S_structdata.th,w,fftw_mag,S_fig);
    elseif strcmp(S_fig.type,'2DfieldFFT') || strcmp(S_fig.type,'2DcoeffsFFT')
        f = plotting_plot2Dfunc(w,fftw_mag,S_fig);
        if S_fig.isdec
            cavity_plotEigenvalues(S_structdata.CFvals,1);
            cavity_plotGainCurve(S_structdata.k_a,S_structdata.g_per);
        end
    end                
end