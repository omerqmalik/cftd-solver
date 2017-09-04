function core_runMicroTDSS(cav_dir,num,pstep,x0,t0,t1,cratio)
    issave_e = [1,1,1,1];
    issave_p = [1,1,1,1];
    issave_d = [0,0,0,0];

    clearvars -global
    addpath(genpath('/tigress/omalik/Time Dynamics/cftd-solver/modules'));
    
    S_coredata = core_init(cav_dir,num,'micro',x0,t0,t1,cratio,pstep);
    
    time_all = zeros(1,pstep);
    tstart_all = tic;
    fprintf('pstep %g\nD0=%f\n',pstep,S_coredata.pump(pstep));
    S_coredata.calc_times = core_calcCoeffs(S_coredata,pstep,issave_e,issave_p,issave_d);
    time_all(pstep)=toc(tstart_all);
    fprintf('Total time: %f\n',time_all);
    
    fprintf('\nSaving all data...');
    
    if issave_e(1) == 1
        rawdata_cleanAll(S_coredata.data_dir,'E','coeffs');
        rawdata_cleanAll(S_coredata.data_dir,'E','field');
    end
    
    if issave_p(1) == 1
        rawdata_cleanAll(S_coredata.data_dir,'P','coeffs');
        rawdata_cleanAll(S_coredata.data_dir,'P','field');
    end
    
    if issave_d(1) == 1
        rawdata_cleanAll(S_coredata.data_dir,'D','coeffs');
    end
    fprintf(' complete.\n');
    
    fprintf('Making figures...');
    userplot_saveMicroFigures(cav_dir,S_coredata.data_dir,num,S_coredata.results_dir,time_all,pstep,S_coredata);
    fprintf(' complete.\n');
end