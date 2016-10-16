function S_coredata = core_init(cav_dir,num,calc_type,varargin)
    load([cav_dir '/TD_init_parameters.mat'],'S_setupdata','S_pumpdata');
    load([cav_dir '/' S_setupdata(num).integral_loc]);
    calc_dir = [cav_dir '/' S_setupdata(num).calc_dir]
    
    %Parameters that go into TDS Solvers
    S_coredata.k_a          = S_setupdata(num).k_a;          %labeled 'ka' elsewhere
    S_coredata.n            = S_setupdata(num).n;            %refractive index
    S_coredata.g_per        = S_setupdata(num).g_per;        %polarization dephasing
    S_coredata.g_par        = S_setupdata(num).g_par;        %inversion relaxation
    S_coredata.CFvals       = S_setupdata(num).CFvals;       %cavity eigenvalues
    S_coredata.A            = A;
    S_coredata.B            = B;
    S_coredata.D0_vec       = D0_vec;
    
    %Parameters required for MATLAB's ODE solver
    S_coredata.tvec         = S_pumpdata.tvec;
    S_coredata.noise_vec    = core_getNoiseVector(S_setupdata(num).nCF,S_setupdata(num).eps);          %magnitude of noise
    
    %Parameters for pump/saving
    S_coredata.basis_type   = S_setupdata(num).basis_type;         %Basis type
    S_coredata.calc_dir     = calc_dir;
    S_coredata.data_dir     = [cav_dir '/' S_setupdata(num).data_dir];
    S_coredata.times_dir    = [cav_dir '/' S_setupdata(num).times_dir];
    S_coredata.results_dir  = [cav_dir '/' S_setupdata(num).results_dir];
    S_coredata.pump_sz      = S_pumpdata.sz;
    S_coredata.sratio       = S_pumpdata.sratio;
    S_coredata.pump         = S_setupdata(num).pump;
    S_coredata.pump_th      = S_setupdata(num).th;
    
    if strcmp(calc_type,'macro')
        S_coredata.calc_times = benchmark_loadTimeFile(S_coredata.times_dir,varargin{1});
        [S_coredata.pump_ind,S_coredata.pumpgrp_ind] = pump_getPumpForTDSSolver(S_coredata.calc_times,S_setupdata(num).pump,S_pumpdata.pg_len,varargin{1});
    end
    
    %Other parameters
    S_coredata.nCF          = S_setupdata(num).nCF;          %# of CF states
    [~,S_coredata.CFvecs]   = cavity_loadBasis([cav_dir '/' S_setupdata(num).basis_loc],S_setupdata(num).k_a);
end

%NOTES (10/9/2016)
%'load' functions should be calling setup_load type functions
%'S_pumpdata' should be phased out
%S_coredata has high similarity to S_setupdata. Reduce to just one?