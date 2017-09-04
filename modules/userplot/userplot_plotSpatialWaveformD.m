function userplot_plotSpatialWaveformD(pstep,pump,t0,x,CFvecs)    
    load('Efield_t.mat');
    [~,tind] = helpers_getClosestMatch(t0,t);
    
    load(['Efield_' num2str(pstep) '.mat']);
    Ecoeffs = Yint(tind);
    
    load(['Pfield_' num2str(pstep) '.mat']);
    Pcoeffs = Yint(tind);
    
    Efield_x = userdata_calcSpatialField(Ecoeffs,CFvecs);
    Pfield_x = userdata_calcSpatialField(Pcoeffs,CFvecs);
    Dfield_x = pump - imag(Efield_x.*conj(Pfield_x));
    
    Efield_x = abs(Efield_x)/max(abs(Efield_x));
    Pfield_x = abs(Pfield_x)/max(abs(Pfield_x));
    Dfield_x = Dfield_x/max(Dfield_x);
    
    plot(x,Efield_x); hold on;
    plot(x,Dfield_x);
end