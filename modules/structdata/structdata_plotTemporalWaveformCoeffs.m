function f = structdata_plotTemporalWaveformCoeffs(S_structdata,S_fig)
    coeffs_t = S_structdata.Y;
    
    coeffs_t = abs(coeffs_t).^2;
    [f,h] = plotting_plot2Dfunc(S_structdata.t,coeffs_t,S_fig);
end