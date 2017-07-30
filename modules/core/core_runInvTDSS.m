function core_runInvTDSS(cav_dir,num,pstep,x0,t1,cratio)
    clearvars -global
    addpath(genpath('/tigress/omalik/Time Dynamics/cftd-solver/modules'));
    
    S_coredata = core_init(cav_dir,num,'inv',x0,0,t1,cratio,pstep);
    
    fprintf('Running calcCoeffs...\n');
    fprintf('pstep %g\nD0=%f\n',pstep,S_coredata.pump(pstep));
    S_coredata.calc_times = core_calcCoeffs(S_coredata,pstep,[1,0,0,1],[1,0,0,1],[0,0,0,0]);
    
    fprintf('Saving CFTD calculation data...');
    rawdata_cleanAll(S_coredata.data_dir,'E','field');
    rawdata_cleanAll(S_coredata.data_dir,'P','field');
    fprintf(' complete.\n\n');
    
    S_EfieldData = structdata_loadMicro(cav_dir,S_coredata.data_dir,num,'E','field',S_coredata.calc_times,pstep,S_coredata.x0);
    S_PfieldData = structdata_loadMicro(cav_dir,S_coredata.data_dir,num,'P','field',S_coredata.calc_times,pstep,S_coredata.x0);
    
    fprintf('Running inversion solver...\n');
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);   %MATLAB defaults: RelTol=1e-3, AbsTol=1e-6 (usually good and fast enough)
    solver_func = @(t,y) MBSolvers_inversion(t,y,S_coredata.pump(pstep),S_coredata.g_par,S_EfieldData.Y,S_PfieldData.Y,S_EfieldData.t);
    
    tref = tic;
    [T,Y] = ode45(solver_func,[0 t1],0,opts);
    inv_time = toc(tref);
    fprintf('Total time: %f\n',inv_time);
    
    fprintf('\nSaving inversion data...');
    rawdata_save(T,Y,S_coredata.data_dir,'D','field',pstep);
    rawdata_cleanAll(S_coredata.data_dir,'D','field');
    fprintf(' complete.\n');
    
    fprintf('Saving plots...');
    userplot_saveEPDFigures(cav_dir,S_coredata.data_dir,num,S_coredata.results_dir,S_coredata.calc_times,inv_time,pstep,S_coredata);
    fprintf(' complete.\n');
end