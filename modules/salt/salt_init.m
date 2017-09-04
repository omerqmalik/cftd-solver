function [S_coredata,S_saltparameters] = salt_init(cav_dir,num,x0,prel_f,sz)
    load([cav_dir '/TD_init_parameters.mat'],'S_setupdata','S_pumpdata');
    load([cav_dir '/' S_setupdata(num).integral_loc]);
    calc_dir = salt_getCalcDir(cav_dir,num,prel_f,sz,S_setupdata);
    
    %Parameters that go into TDS Solvers
    S_coredata.k_a          = S_setupdata(num).k_a;          %labeled 'ka' elsewhere
    S_coredata.n            = S_setupdata(num).n;            %refractive index
    S_coredata.g_per        = S_setupdata(num).g_per;        %polarization dephasing
    S_coredata.g_par        = S_setupdata(num).g_par;        %inversion relaxation
    S_coredata.CFvals       = S_setupdata(num).CFvals;       %cavity eigenvalues
    
    %Parameters for pump/saving
    S_coredata.basis_type   = S_setupdata(num).basis_type;         %Basis type
    S_coredata.basis_loc    = S_setupdata(num).basis_loc;
    S_coredata.calc_dir     = calc_dir;
    S_coredata.pump_th      = S_setupdata(num).th;
    
    S_coredata.prel_f  = prel_f;
    S_coredata.pump_sz = sz;
    
    %Other parameters
    S_coredata.nCF          = S_setupdata(num).nCF;          %# of CF states
    [~,S_coredata.CFvecs]   = cavity_loadBasis([cav_dir '/' S_setupdata(num).basis_loc],S_setupdata(num).k_a);
    S_coredata.x0           = x0;
    
    %Set up directories
    [~,data_dir] = fileparts(S_setupdata(num).data_dir);
    [~,results_dir] = fileparts(S_setupdata(num).results_dir);
    
    S_coredata.data_dir     = [calc_dir '/' data_dir];
    S_coredata.results_dir  = [calc_dir '/' results_dir];
%         S_coredata.cp_dir       = [pre_dir '/' S_setupdata(num).cp_dir];
%         S_coredata.diag_dir     = [pre_dir '/' S_setupdata(num).diag_dir];
%         S_coredata.temp_dir     = [getenv('CFTD_TEMP_PATH') '/' S_setupdata(num).temp_dir];
%         S_coredata.times_dir    = [pre_dir '/' S_setupdata(num).times_dir];
    
    mkdir(S_coredata.data_dir);
    mkdir(S_coredata.results_dir);
        
    SALT_parameters; 
    
    S_saltparameters.type   = SALT_type;
    
    S_saltparameters.skCtr  = S_coredata.k_a;               %Central k-value around which 'slices' are taken for calculation
    S_saltparameters.skWnd  = TN_skWnd*S_coredata.g_per;    %Window of k-values within which to take slices
    S_saltparameters.skN    = TN_skN;                       %Number of slices to take
    S_saltparameters.sskN   = TN_sskN;                      %Number of subslices to take
    
    S_saltparameters.ntest  = TN_ntest;                     %An arbitrary large number of iterations for fine tuning k and threshold
    S_saltparameters.kerr   = TN_kerr;                      %Cease iterations When iterative corrections to k become smaller than this value
    S_saltparameters.therr  = TN_therr;                     %Cease iterations when iterative corrections to threshold become smaller than this value
    
    S_saltparameters.fskWnd = TN_fskWnd;                    %Window of k-values within which to take fine slices (.006)
    S_saltparameters.fskN   = TN_fskN;                      %Number of fine slices to take, previously 200 (100)
    S_saltparameters.kcri   = TN_kcri;                      %Only when the change in k after recalculating the B matrix is smaller than kcri, we do fine tuning
    S_saltparameters.kprec  = TN_kprec;                     %If the change in k by reducing the step length('stepdownF') is smaller than dkprec,
    S_saltparameters.dkprec = TN_dkprec;                    %or if the stepdownF is smaller than kprec, precision goal is achieved
    
    %mod_TSfsolve
    S_saltparameters.dtves      = TN_dtves;        %small differential number to multiply with tves
    S_saltparameters.upperbound = TN_upperbound;   %If the new-found direction vector's largest element is smaller than this value, then normalize

    %mod_MultimodeFsolve
    S_saltparameters.MMfslvOptTolFun = TN_MMfslvOptTolFun;   %Tolerance on the function value returned by fsolve
    S_saltparameters.MMfslvOptTolX   = TN_MMfslvOptTolX;     %Tolerance on the iterative seed x
    
    save([calc_dir '/setup_parameters.mat'],'S_coredata','S_saltparameters');
end