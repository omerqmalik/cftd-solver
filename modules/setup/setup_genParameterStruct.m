function [S_setupdata,S_pumpdata] = setup_genParameterStruct(param_vecs,S_pumpdata,basis_type,rframe)
    t0        = S_pumpdata.t0;
    t1        = S_pumpdata.t1;
    cratio    = S_pumpdata.cratio;
    prel_i    = S_pumpdata.prel_i;
    prel_f    = S_pumpdata.prel_f;
    p_dens    = S_pumpdata.p_dens;
    pg_len    = S_pumpdata.pg_len;
    xdens     = S_pumpdata.xdens;
    
    if isfield(S_pumpdata,'type')
        pump_type = S_pumpdata.type;
    else
        pump_type = 'simple';
        S_pumpdata.type = pump_type;
    end
    
    %Pump details
    if prel_f == prel_i
        S_pumpdata.sz = 1;
        pumprel = prel_i;
    else
        S_pumpdata.sz = ceil((prel_f - prel_i)*p_dens);
        pumprel = linspace(prel_i,prel_f,S_pumpdata.sz);
        if strcmp(pump_type,'hysteresis')
            pumprel = [pumprel fliplr(pumprel(1:(end-1)))];
            S_pumpdata.sz = length(pumprel);
            pg_len = S_pumpdata.sz;
        end
    end
    S_pumpdata.pgroups = ceil(S_pumpdata.sz/pg_len);
    pgroups = S_pumpdata.pgroups;
    sz = S_pumpdata.sz;
    
    %Time vector details
    time_vec = setup_divideTimeVector(t0,t1,cratio);
    S_pumpdata.tvec = time_vec;
    
    % Initialize empty structs
    leN = size(param_vecs,1);
    S_setupdata = struct('k_a',[],'n',[],'nCF',[],'g_per',[], ...
        'g_par',[],'prel_i',[],'prel_f',[],'pump',[],'th',[], ...
        'eps',[],'dx',[],'x',[],'w_FSR',[],'Na',[],'M',[],'rframe',[], ...
        'CFvals',[],'basis_type',[],'basis_loc',[],'integral_loc',[], ...
        'calc_path',[],'calc_dir',[],'data_dir',[],'times_dir',[], ...
        'results_dir',[]);
    S_setupdata(leN).nCF = [];
    
    integral_dir = 'integrals';
    mkdir(integral_dir);

    index = 1;
    done_NCF = [];
    done_k_a = [];
    done_n   = [];
    for i = 1:leN
        k_a   = param_vecs(i,1);
        n     = param_vecs(i,2);
        nCF   = param_vecs(i,3);
        g_per = param_vecs(i,4);
        g_par = param_vecs(i,5);
        eps   = param_vecs(i,6);
        
        if nargin < 4
            rframe = k_a;
        end
        
        %Name directories
        calc_dir    = ['num_' num2str(index)];
        data_dir    = [calc_dir '/data'];
        times_dir   = [calc_dir '/times'];
        results_dir = [calc_dir '/results'];
        basis_dir   = [calc_dir '/bases'];
        cp_dir      = [calc_dir '/checkpoints'];
        diag_dir    = [calc_dir '/diagnostics'];
        
        [~,tempdirname] = fileparts(tempname);
        temp_dir = tempdirname;
        
        %Create directories
        mkdir(data_dir);
        mkdir(times_dir);
        mkdir(results_dir);
        mkdir(basis_dir);
        mkdir(cp_dir);
        mkdir(diag_dir);
        mkdir([getenv('CFTD_TEMP_PATH') '/' temp_dir]);
        
        %Calculate/get basis
        [CFvals,CFvecs,dx,nx,x,w_FSR,Na,M] = cavity_calcModes(basis_type,k_a,n,nCF,xdens);
        basis_loc = [basis_dir '/' basis_type 'basis_nCF_' num2str(nCF) '_ka_' num2str(k_a) '_n_' num2str(n) '.mat'];
        save(basis_loc,'CFvals','CFvecs','basis_type','dx','nx','x','w_FSR','Na','M');

        %Calculate integrals (only once)
        if sum((done_NCF == nCF) & (done_k_a == k_a) & (done_n == n)) == 0
            integral_loc = [integral_dir '/TD_integrals_nCF_' num2str(nCF) '_ka_' num2str(k_a) '_n_' num2str(n) '.mat'];
            if ~exist(integral_loc,'file')
                int_time = tic;
                [A,A_index,B,D0_vec,D_index] = setup_calcOverlapIntegrals(nCF,n^2,ones(length(x),1),CFvecs,basis_type,dx,M);
                int_time = toc(int_time);
                fprintf(['nCF = ' num2str(nCF) ': %fs\n'],int_time);

                save(integral_loc,'A','A_index','B','D0_vec','D_index','-v7.3');
            end
    
            done_NCF = [done_NCF nCF];
            done_k_a = [done_k_a k_a];
            done_n   = [done_n n];
        else
            integral_loc = [integral_dir '/TD_integrals_nCF_' num2str(nCF) '_ka_' num2str(k_a) '_n_' num2str(n) '.mat'];
        end

        %Get threshold
        nu0 = real(CFvals(1));
        kappa0 = -imag(CFvals(1));
        th = cavity_calcDth(n,g_per,nu0,kappa0,k_a);
        
        %Create pump vector
        pump = th*pumprel;

        %Populate struct
        S_setupdata(index).k_a          = k_a;
        S_setupdata(index).n            = n;
        S_setupdata(index).nCF          = nCF;
        S_setupdata(index).g_per        = g_per;
        S_setupdata(index).g_par        = g_par;
        S_setupdata(index).eps          = eps;
        S_setupdata(index).prel_i       = prel_i;
        S_setupdata(index).prel_f       = prel_f;
        S_setupdata(index).pump         = pump;
        S_setupdata(index).th           = th;
        S_setupdata(index).dx           = dx;
        S_setupdata(index).x            = x;
        S_setupdata(index).w_FSR        = w_FSR;
        S_setupdata(index).Na           = Na;
        S_setupdata(index).M            = M;
        S_setupdata(index).rframe       = rframe;
        S_setupdata(index).CFvals       = CFvals;
        S_setupdata(index).basis_type   = basis_type;
        S_setupdata(index).calc_path    = pwd;
        S_setupdata(index).calc_dir     = calc_dir;
        S_setupdata(index).data_dir     = data_dir;
        S_setupdata(index).times_dir    = times_dir;
        S_setupdata(index).results_dir  = results_dir;
        S_setupdata(index).temp_dir     = temp_dir;
        S_setupdata(index).cp_dir       = cp_dir;
        S_setupdata(index).diag_dir     = diag_dir;
        S_setupdata(index).basis_loc    = basis_loc;
        S_setupdata(index).integral_loc = integral_loc;
        
        %Initialize noise vectors
        noise_vec = setup_getNoiseVector(nCF,eps);
        for ii = 1:length(pumprel)
            checkpoints_fn = core_getCheckpointFn(ii,cp_dir);
            core_saveCheckpoints(time_vec(1),noise_vec,checkpoints_fn);
            if strcmp(pump_type,'hysteresis')
                break;
            end
        end
        
        % Create empty files for storing calculation times
        benchmark_createTimeFiles(times_dir,pgroups,pg_len,sz,length(time_vec));

        index = index + 1;
    end
end