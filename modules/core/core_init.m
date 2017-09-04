function S_coredata = core_init(cav_dir,num,calc_type,x0,varargin)
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
    
    %Parameters for pump/saving
    S_coredata.basis_type   = S_setupdata(num).basis_type;         %Basis type
    S_coredata.basis_loc    = S_setupdata(num).basis_loc;
    S_coredata.calc_dir     = calc_dir;
    S_coredata.pump_th      = S_setupdata(num).th;
    S_coredata.pump         = S_setupdata(num).pump;
    S_coredata.pump_type    = S_pumpdata.type;
    
    %Other parameters
    S_coredata.nCF          = S_setupdata(num).nCF;          %# of CF states
    [~,S_coredata.CFvecs]   = cavity_loadBasis([cav_dir '/' S_setupdata(num).basis_loc],S_setupdata(num).k_a);
    S_coredata.x0           = x0;
    
    %Parameters used by macro calculation
    if strcmp(calc_type,'macro')
        pgroup = varargin{1};
        pre_dir = cav_dir;
        
        %Parameters required for MATLAB's ODE solver
        S_coredata.tvec         = S_pumpdata.tvec;
        
        %Set up directories
        S_coredata.data_dir     = [pre_dir '/' S_setupdata(num).data_dir];
        S_coredata.times_dir    = [pre_dir '/' S_setupdata(num).times_dir];
        S_coredata.results_dir  = [pre_dir '/' S_setupdata(num).results_dir];
        S_coredata.cp_dir       = [pre_dir '/' S_setupdata(num).cp_dir];
        S_coredata.diag_dir     = [pre_dir '/' S_setupdata(num).diag_dir];
        S_coredata.temp_dir     = [getenv('CFTD_TEMP_PATH') '/' S_setupdata(num).temp_dir];
        
        %Set up pump
        S_coredata.pump_sz      = S_pumpdata.sz;
        S_coredata.sratio       = S_pumpdata.sratio;
        
        %Get starting pump step/index and pumpgroup
        S_coredata.calc_times = benchmark_loadTimeFile(S_coredata.times_dir,pgroup);
        [S_coredata.pump_ind,S_coredata.pumpgrp_ind] = pump_getPumpForTDSSolver(S_coredata.calc_times,S_setupdata(num).pump,S_pumpdata.pg_len,S_pumpdata.type,varargin{1});
    elseif strcmp(calc_type,'micro') || strcmp(calc_type,'inv')
        t0     = varargin{1};
        t1     = varargin{2};
        cratio = varargin{3};
        pstep  = varargin{4};
        
        pre_dir  = [calc_dir '/' calc_type '/p_' num2str(pstep) '_t0_' num2str(t0) '_t1_' num2str(t1)];
        
        [~,data_dir]    = fileparts(S_setupdata(num).data_dir);
        [~,times_dir]   = fileparts(S_setupdata(num).times_dir);
        [~,results_dir] = fileparts(S_setupdata(num).results_dir);
        [~,cp_dir]      = fileparts(S_setupdata(num).cp_dir);
        [~,diag_dir]    = fileparts(S_setupdata(num).diag_dir);
        [~,temp_dir]    = fileparts(tempname);
        
        data_dir        = [pre_dir '/' data_dir];
        times_dir       = [pre_dir '/' times_dir];
        results_dir     = [pre_dir '/' results_dir];
        cp_dir          = [pre_dir '/' cp_dir];
        diag_dir        = [pre_dir '/' diag_dir];
        temp_dir        = [getenv('CFTD_TEMP_PATH') '/' temp_dir];
        
        mkdir(pre_dir);
        mkdir(times_dir);
        mkdir(data_dir);
        mkdir(results_dir);
        mkdir(cp_dir);
        mkdir(diag_dir);
        mkdir(temp_dir);
        
        S_coredata.data_dir     = data_dir;
        S_coredata.times_dir    = times_dir;
        S_coredata.results_dir  = results_dir;
        S_coredata.cp_dir       = cp_dir;
        S_coredata.diag_dir     = diag_dir;
        S_coredata.temp_dir     = temp_dir;
        
        S_coredata.pump_sz      = 1;
        S_coredata.sratio       = 1;
        
        [t,noise_vec] = core_loadCheckpoints(core_getCheckpointFn(pstep,[cav_dir '/' S_setupdata(num).cp_dir]));
        S_coredata.tvec = setup_divideTimeVector(t(t == max(t(t<=t0))),t1,cratio);
        core_saveCheckpoints(S_coredata.tvec(1),noise_vec(:,t == max(t(t<=t0))),core_getCheckpointFn(pstep,cp_dir));
        
        diag_ind = (1:(S_coredata.nCF+1):S_coredata.nCF^2)+2*S_coredata.nCF;
        
        %To select random indices
%         all_ind  = (2*S_coredata.nCF+1):(S_coredata.nCF^2+2*S_coredata.nCF);
%         diff_ind = setdiff(all_ind,diag_ind);
%         rand_ind = diff_ind(randperm(length(diff_ind)));
%         S_coredata.Dsave = [diag_ind rand_ind(1:S_coredata.nCF)];
        
        %To select the +-1 diagonal indices
        diagp_ind = diag_ind(2:end)-1;
        diagm_ind = diag_ind(1:end-1)+1;
        S_coredata.Dsave = [diag_ind diagp_ind diagm_ind];
    end
end

%NOTES (10/9/2016)
%'load' functions should be calling setup_load type functions
%'S_pumpdata' should be phased out
%S_coredata has high similarity to S_setupdata. Reduce to just one?