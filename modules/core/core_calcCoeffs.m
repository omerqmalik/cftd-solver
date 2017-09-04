function calc_times = core_calcCoeffs(S_coredata,pstep,issave_e,issave_p,issave_d)
    global benchmarking;
    
    fID = fopen('mem_file.txt','a');
    
    
    
    %Get first group from tvec to begin saving
    tvec        = S_coredata.tvec;
    len_tvec    = length(tvec);

    %Initialize
    [t_last,Y_last] = core_loadCheckpoints(core_getCheckpointFn(pstep,S_coredata.cp_dir));
    noise_vec  = Y_last(:,end);
    t_ind = find(tvec == t_last(end));
    
    if (benchmarking) 
        if t_ind == 1
            calc_times = zeros(len_tvec-1,1);
            benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
        else
            calc_times = benchmark_loadTimeForPstep(S_coredata.times_dir,pstep);
        end
    end
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);   %MATLAB defaults: RelTol=1e-3, AbsTol=1e-6 (usually good and fast enough)
    
    %Run ODE solver, save E_t, and store T and Y in memory to return at end
    for j = t_ind:(len_tvec - 1)
        tstart_in = tic;
        [T,Y] = ode45Wrapper(S_coredata, pstep, [tvec(j) tvec(j+1)],noise_vec,opts);
        calc_times(j) = toc(tstart_in);
        noise_vec = Y(end,:);
        
        rawdata_save([T(1) T(end)],mean(abs(Y(:,(2*S_coredata.nCF+1):end))).',S_coredata.temp_dir,'D','avgabs',pstep,j);
        if issave_d(1) == 1
            rawdata_save(T,Y(:,S_coredata.Dsave),S_coredata.temp_dir,'D','coeffs',pstep,j);
        end
        Y = Y(:,1:2*S_coredata.nCF);

        %Save e_m and p_m
        if issave_e(1) == 1
            rawdata_save(T,Y(:,1:S_coredata.nCF),S_coredata.temp_dir,'E','coeffs',pstep,j);
        end
        
        if issave_p(1) == 1
            rawdata_save(T,Y(:,(S_coredata.nCF+1):2*S_coredata.nCF),S_coredata.temp_dir,'P','coeffs',pstep,j);
        end
        if (benchmarking)        
            %bookkeeping
            fprintf('iteration %g: %fs, memory: %d\n',j,calc_times(j),java.lang.Runtime.getRuntime.totalMemory);
            benchmark_saveTimeForPstep(S_coredata.times_dir,pstep,calc_times);
            core_saveCheckpoints(tvec(j+1),noise_vec.',core_getCheckpointFn(pstep,S_coredata.cp_dir));
        end
    end
    if (benchmarking) 
        fprintf('Total time: %fs\n\n',sum(calc_times));
    end
    
    if issave_e(1) == 1
        diag_save1DCalcData(S_coredata,'E','coeffs',pstep,issave_e(2),issave_e(3),issave_e(4));
    end
    
    if issave_p(1) == 1
        diag_save1DCalcData(S_coredata,'P','coeffs',pstep,issave_p(2),issave_p(3),issave_p(4));
    end
    
    if issave_d(1) == 1
        diag_save1DCalcData(S_coredata,'D','coeffs',pstep,issave_d(2),issave_p(3),issave_p(4));
        diag_save2DCalcData(S_coredata,'D','avgabs',pstep);
    end

    fclose(fID);
end