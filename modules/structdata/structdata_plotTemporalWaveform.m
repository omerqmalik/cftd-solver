function f = structdata_plotTemporalWaveform(S_structdata,S_fig)
    if strcmp(S_structdata.type,'coeffs')
        [~,CFvecs,~,x] = cavity_loadBasis(S_structdata.basis_loc,S_structdata.k_a);
        field_t    = userdata_calcTemporalField(S_structdata.Y,CFvecs,x);
    elseif strcmp(S_structdata.type,'field')
        field_t = S_structdata.Y;
    end
    
    field_t = abs(field_t);
    [f,h] = plotting_plot2Dfunc(S_structdata.t,field_t,S_fig);
    if length(h) > 1
        c = jet(length(h));
        for i = 1:length(h)
            h(i).Color = c(i,:);
        end
    end
end