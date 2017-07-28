function core_runTDSS(cav_dir,num,pgroup)
    clearvars -global
    global benchmarking;
    global usingc;
    global save_ode_params;
    usingc = true;
    benchmarking = true;
    save_ode_params = true;
    addpath(genpath('/tigress/omalik/Time Dynamics/cftd-solver/modules'));
    
    fprintf('num: %d\npgroup: %d\n',num,pgroup);
    S_coredata = core_init(cav_dir,num,'macro',pgroup);
    set(0,'DefaultFigureVisible','off');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 TIME DYNAMICAL CALCULATION STARTS HERE                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (benchmarking)
        tstart_all = tic;
    end
    for i = 1:length(S_coredata.pump_ind)
        %initialize
        pstep  = S_coredata.pump_ind(i);
        pgstep = S_coredata.pumpgrp_ind(i);
        if (benchmarking) 
            fprintf('pstep %g\nD0=%f\n',pstep,S_coredata.pump(pstep));
        end
        
        %calculate and save
        [S_coredata.calc_times(:,pgstep)] = core_calcCoeffs(S_coredata,pstep,1,1,0);
        
        if strcmp(S_coredata.pump_type,'hysteresis') && i < length(S_coredata.pump_ind)
            [~,Y_last] = core_loadCheckpoints(core_getCheckpointFn(pstep,S_coredata.cp_dir));
            core_saveCheckpoints(S_coredata.tvec(1),Y_last(:,end),core_getCheckpointFn(pstep+1,S_coredata.cp_dir));
        end
        if (benchmarking)
            %save benchmark
            benchmark_saveTimeFile(S_coredata.times_dir,pgroup,S_coredata.calc_times);
            benchmark_saveDoneFile(S_coredata.times_dir,pstep);
        end
    end
    if (benchmarking)
        time_all=toc(tstart_all);
        fprintf('Total time: %f\n',time_all);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  TIME DYNAMICAL CALCULATION ENDS HERE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Post processing
    datafiles = dir([S_coredata.times_dir '/done_*.mat']);
    
    if (length(datafiles) == S_coredata.pump_sz)
        fprintf('\nSaving all data...');
        rawdata_cleanAll(S_coredata.data_dir,'E','coeffs');
        rawdata_cleanAll(S_coredata.data_dir,'E','field');
        %rawdata_cleanAll(S_coredata.data_dir,'P','field');
        benchmark_consolidateTimes(S_coredata.times_dir);
        fprintf(' complete.\n');
        
        fprintf('Making figures...');
        userplot_saveFigures(cav_dir,num,S_coredata.results_dir);
        fprintf(' complete.\n');
    end
    set(0,'DefaultFigureVisible','on');
end

%Things to plot/save:
%1. Field 3Dfft
%2. Field 2Dfft at each power (optional?)
%3. Field 2Dfft movie
%4. Field time signal at each power (abs only)
%5. Field time signal movie (abs only)
