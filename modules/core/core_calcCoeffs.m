function [T,Y,calc_times] = core_calcCoeffs(S_coredata,pstep)

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
    end
    
    %Get first group from tvec to begin saving
    tvec        = S_coredata.tvec;
    sratio      = S_coredata.sratio;
    len_tvec    = length(tvec);
    save_group0 = get_save_group0(tvec,sratio);

    %Initialize
    noise_vec  = S_coredata.noise_vec;
    calc_times = zeros(len_tvec-1,1);
    benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
    T_temp = [];
    Y_temp = [];
    is_saving = 0;
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);   %MATLAB defaults: RelTol=1e-3, AbsTol=1e-6 (usually good enough)
    
    %Run ODE solver, save E_t, and store T and Y in memory to return at end
    for j = 1:(len_tvec - 1)
        tstart_in = tic;
        [T,Y] = ode45(solver_func, [tvec(j) tvec(j+1)],noise_vec,opts);
        calc_times(j) = toc(tstart_in);
        noise_vec = Y(end,:);

        %Save E_t
        E_t_next = userdata_calcTemporalField(Y(:,1:S_coredata.nCF),S_coredata.CFvecs);
        T_next   = T;
        if exist(rawdata_getFileName(S_coredata.data_dir,'E','field',pstep),'file') == 2
            [T,E_t] = rawdata_load(S_coredata.data_dir,'E','field',pstep);
            E_t = [E_t; E_t_next];
            T   = [T; T_next];
            rawdata_save(T,E_t,S_coredata.data_dir,'E','field',pstep);
        else
            E_t = E_t_next;
            rawdata_save(T,E_t,S_coredata.data_dir,'E','field',pstep);
        end
        clear E_t_next;
        
        %Save P_t
        P_t_next = userdata_calcTemporalField(Y(:,(S_coredata.nCF+1):2*S_coredata.nCF),S_coredata.CFvecs);
        if exist(rawdata_getFileName(S_coredata.data_dir,'P','field',pstep),'file') == 2
            [~,P_t] = rawdata_load(S_coredata.data_dir,'P','field',pstep);
            P_t = [P_t; P_t_next];
            rawdata_save(T,P_t,S_coredata.data_dir,'P','field',pstep);
        else
            P_t = P_t_next;
            rawdata_save(T,P_t,S_coredata.data_dir,'P','field',pstep);
        end
        clear P_t_next;
        
        %Appropriately store T and Y in memory
        if j < save_group0
            clear T Y T_next E_t P_t;
        else
            if is_saving == 0
                Y      = Y(T_next >= tvec(end)*(1-sratio),:);
                T_next = T_next(T_next >= tvec(end)*(1-sratio));
                is_saving = 1;
            end
            T_temp = [T_temp; T_next(2:end)];
            Y_temp = [Y_temp; Y(2:end,:)];
            if j < (len_tvec - 1)
                clear T_next E_t P_t;
            end
        end
        fprintf('iteration %g: %fs\n',j,calc_times(j));
        benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
    end
    fprintf('Total time: %fs\n\n',sum(calc_times));     %Is the data below being saved 'clean'?
    
    T = T_temp;
    Y = Y_temp;
end

function save_group0 = get_save_group0(tvec,sratio)
    save_t0 = round(tvec(end)*(1 - sratio));
    tvec_new = sort([tvec save_t0]);
    pos_new = find(tvec_new == save_t0);
    if length(pos_new) > 1
        save_group0 = pos_new(1);
    else
        if pos_new <= 2
            save_group0 = 1;
        elseif pos_new >= length(tvec)
            save_group0 = length(tvec) - 1;
        else
            save_group0 = pos_new - 1;
        end
    end
end
