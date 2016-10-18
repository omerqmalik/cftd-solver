function core_runTDSS(cav_dir,num,pgroup)
    clearvars -global
    addpath(genpath('/tigress/omalik/Time Dynamics/cftd-solver/modules'));
    
    fprintf('num: %d\npgroup: %d\n',num,pgroup);
    S_coredata = core_init(cav_dir,num,'macro',pgroup);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 TIME DYNAMICAL CALCULATION STARTS HERE                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tstart_all = tic;
    for i = 1:length(S_coredata.pump_ind)
        %initialize
        pstep  = S_coredata.pump_ind(i);
        pgstep = S_coredata.pumpgrp_ind(i);
        fprintf('pstep %g\nD0=%f\n',pstep,S_coredata.pump(pstep));
        
        %calculate
        [T,Y,S_coredata.calc_times(:,pgstep)] = core_calcCoeffs(S_coredata,pstep);
        
        %save
        Y = Y(:,1:S_coredata.nCF);
        rawdata_save(T,Y,S_coredata.data_dir,'E','coeffs',pstep)
        benchmark_saveTimeFile(S_coredata.times_dir,pgroup,S_coredata.calc_times);
    end
    time_all=toc(tstart_all);
    fprintf('Total time: %f\n',time_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  TIME DYNAMICAL CALCULATION ENDS HERE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Post processing
    datafiles = dir([S_coredata.data_dir '/Ecoeffs_*.mat']);
    
    if length(datafiles) == S_coredata.pump_sz
        fprintf('\nSaving all data...');
        rawdata_cleanAll(S_coredata.data_dir,'E','coeffs');
        rawdata_cleanAll(S_coredata.data_dir,'E','field');
        rawdata_cleanAll(S_coredata.data_dir,'P','field');
        benchmark_consolidateTimes(S_coredata.times_dir);
        fprintf(' complete.\n');
        
        fprintf('Making figures...');
        structdata_saveFigures(cav_dir,num,S_coredata.results_dir);
        fprintf(' complete.\n');
    end
end

%Things to plot/save:
%1. Field 3Dfft
%2. Field 2Dfft at each power (optional?)
%3. Field 2Dfft movie
%4. Field time signal at each power (abs only)
%5. Field time signal movie (abs only)
