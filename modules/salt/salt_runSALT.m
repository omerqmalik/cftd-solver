function salt_runSALT(cav_dir,num,prel_f,sz)
    clearvars -global
    
    x0 = 0.327;
    
    [S_coredata,S_saltparameters] = salt_init(cav_dir,num,x0,prel_f,sz);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SALT CALCULATION STARTS HERE                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S_rough = salt_scanRoughThresholds(S_coredata,S_saltparameters);
    save([S_coredata.data_dir '/rough_scan.mat'],'S_rough');
    
    S_fine  = salt_refineThresholds(S_rough,S_coredata,S_saltparameters,zeros(size(S_coredata.CFvecs,1),1));
    save([S_coredata.data_dir '/fine_scan.mat'],'S_fine');
    
    salt_calcNLSolution(S_fine,S_coredata,S_saltparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  TIME DYNAMICAL CALCULATION ENDS HERE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end