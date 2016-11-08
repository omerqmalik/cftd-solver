function S_structdata = structdata_load(cav_dir,num,data_id,data_type,chunk,varargin)
    S_setupdata = setup_loadParameters(cav_dir,num);
    calc_times  = benchmark_loadTimes([cav_dir '/' S_setupdata.times_dir]);
    
    if isempty(varargin)
        psteps = 1:length(S_setupdata.pump);
    else
        pump0 = S_setupdata.th*varargin{1};
        if length(varargin) == 1
            pump1 = pump0;
        elseif length(varargin) == 2
            pump1 = S_setupdata.th*varargin{2};
        end
        psteps = pump_getPumpSteps(S_setupdata.pump,pump0,pump1);
    end
    psteps = psteps(sum(calc_times(:,psteps),1) > 0);
    
    [t,Y]  = userdata_load([cav_dir '/' S_setupdata.data_dir],data_id,data_type,psteps);
    t = helpers_truncVector(t.',chunk);
    Y = helpers_truncVector(Y,chunk);

    S_structdata.t           = t;
    S_structdata.Y           = Y;
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
    
    if strcmp(S_structdata.basis_type,'RING')
        S_structdata.Na    = S_setupdata.Na;
        S_structdata.M     = S_setupdata.M;
        S_structdata.w_FSR = S_setupdata.w_FSR;
    end
end