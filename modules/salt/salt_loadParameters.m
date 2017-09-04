function [S_coredata,S_saltparameters] = salt_loadParameters(cav_dir,num,prel_f,sz)
    load([cav_dir '/TD_init_parameters.mat'],'S_setupdata');
    calc_dir = salt_getCalcDir(cav_dir,num,prel_f,sz,S_setupdata);
    
    load([calc_dir '/setup_parameters.mat']);
end