function S_structdata = structdata_getGeneralProperties(data_id,data_type,S_setupdata,calc_times,psteps)
    S_structdata.id          = data_id;
    S_structdata.type        = data_type;
    S_structdata.calc_times  = calc_times(:,psteps);
    S_structdata.pump        = S_setupdata.pump(psteps);
    S_structdata.k_a         = S_setupdata.k_a;
    S_structdata.n           = S_setupdata.n;
    S_structdata.nCF         = S_setupdata.nCF;
    S_structdata.g_per       = S_setupdata.g_per;
    S_structdata.g_par       = S_setupdata.g_par;
    S_structdata.eps         = S_setupdata.eps;
    S_structdata.th          = S_setupdata.th;
    S_structdata.rframe      = S_setupdata.rframe;
    S_structdata.CFvals      = S_setupdata.CFvals;
    S_structdata.basis_type  = S_setupdata.basis_type;
    S_structdata.basis_loc   = S_setupdata.basis_loc;   
    S_structdata.pump_type   = S_setupdata.pump_type;
    
    if strcmp(S_structdata.basis_type,'RING')
        S_structdata.Na    = S_setupdata.Na;
        S_structdata.M     = S_setupdata.M;
        S_structdata.w_FSR = S_setupdata.w_FSR;
    end
end