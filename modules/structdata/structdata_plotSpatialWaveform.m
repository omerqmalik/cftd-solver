function f = structdata_plotSpatialWaveform(S_structdata,S_fig)
    [~,CFvecs,~,x] = cavity_loadBasis(S_structdata.basis_loc,S_structdata.k_a);
    field_x    = userdata_calcSpatialField(S_structdata.Y,CFvecs);
    
    [f,h] = plotting_plot2Dfunc(x,abs(field_x),S_fig);
%     plotting_plot2Dfunc(x,real(field_x),S_fig);
%     plotting_plot2Dfunc(x,imag(field_x),S_fig);
end