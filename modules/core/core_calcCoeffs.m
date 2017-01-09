function calc_times = core_calcCoeffs(S_coredata,pstep,issave_e,issave_p,issave_d)

    %Generate appropriate function handle
    g_per      = S_coredata.g_per;
    g_par      = S_coredata.g_par;
    k_a        = S_coredata.k_a;
    CFvals     = S_coredata.CFvals;
    A          = S_coredata.A;
    B          = S_coredata.B;
    pump       = S_coredata.pump;
    D0_vec     = S_coredata.D0_vec;
    n          = S_coredata.n;
    basis_type = S_coredata.basis_type;
    
    pump_pwr = pump(pstep)*D0_vec;
    if strcmp(basis_type,'RING')
        solver_func = @(t,y) TDSSolvers_RING(t,y,g_per,g_par,k_a,CFvals,A,pump_pwr,length(CFvals),n);
    elseif strcmp(basis_type,'FP')
        solver_func = @(t,y) TDSSolvers_FP(t,y,g_per,g_par,k_a,CFvals,A,pump_pwr,length(CFvals),n);
    elseif strcmp(basis_type,'UCF')
        solver_func = @(t,y) TDSSolvers_UCF(t,y,g_per,g_par,k_a,CFvals,A,B,pump_pwr,length(CFvals));
    end
    
    %Get first group from tvec to begin saving
    tvec        = S_coredata.tvec;
    len_tvec    = length(tvec);

    %Initialize
    [t_last,Y_last] = core_loadCheckpoints(core_getCheckpointFn(pstep,S_coredata.cp_dir));
    noise_vec  = Y_last(:,end);
    t_ind = find(tvec == t_last(end));
    if t_ind == 1
        calc_times = zeros(len_tvec-1,1);
        benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
    else
        calc_times = benchmark_loadTimeForPstep(S_coredata.times_dir,pstep);
    end
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);   %MATLAB defaults: RelTol=1e-3, AbsTol=1e-6 (usually good and fast enough)
    
    %Run ODE solver, save E_t, and store T and Y in memory to return at end
    for j = t_ind:(len_tvec - 1)
        tstart_in = tic;
        [T,Y] = ode45(solver_func, [tvec(j) tvec(j+1)],noise_vec,opts);
        calc_times(j) = toc(tstart_in);
        noise_vec = Y(end,:);
        
        rawdata_save([T(1) T(end)],mean(abs(Y(:,(2*S_coredata.nCF+1):end))).',S_coredata.temp_dir,'D','avgabs',pstep,j);
        if issave_d
            rawdata_save(T,Y(:,S_coredata.Dsave),S_coredata.temp_dir,'D','coeffs',pstep,j);
        end
        Y = Y(:,1:2*S_coredata.nCF);

        %Save e_m and p_m
        if issave_e
            rawdata_save(T,Y(:,1:S_coredata.nCF),S_coredata.temp_dir,'E','coeffs',pstep,j);
        end
        
        if issave_p
            rawdata_save(T,Y(:,(S_coredata.nCF+1):2*S_coredata.nCF),S_coredata.temp_dir,'P','coeffs',pstep,j);
        end
                
        %bookkeeping
        fprintf('iteration %g: %fs\n',j,calc_times(j));
        benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
        core_saveCheckpoints(tvec(j+1),noise_vec.',core_getCheckpointFn(pstep,S_coredata.cp_dir));
    end
    fprintf('Total time: %fs\n\n',sum(calc_times));
    
    if issave_e
        diag_save1DCalcData(S_coredata,'E','coeffs',pstep,1,1,1);
    end
    
    if issave_p
        diag_save1DCalcData(S_coredata,'P','coeffs',pstep,0,1,0);
    end
    
    if issave_d
        diag_save1DCalcData(S_coredata,'D','coeffs',pstep,1,0,0);
    end
    diag_save2DCalcData(S_coredata,'D','avgabs',pstep);
end