clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SET VECTOR PARAMETERS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Variables prepended A_ are arrays. You may specify more than one
% value for parameters with A_ prepended.

%This file is for preparing  Fabry-Perot cavities for calculations
%Some options below are unique to this basis type
basis_type = 'FP';

%Specify the basis vectors to generate:
%Basis vectors are calculated in a frequency window centered at a value
%defined by cmode (central mode). cmode_type can be either 'Na' or 'ka'.
%When cmode_type is set to 'ka', then A_basis_cmode must specify exact
%value (or values) for the center frequency. If cmode_type is set to 'Na'
%then the index (or indices) of the cavity mode you want to use as the
%center frequency should be specified. The corresponding frequency value
%will be the center of the frequency window and it will also be the value
%of the atomic frequency. See examples:

% %Select mode frequency:
% %Gain center may not lie on top of a cavity mode
% cmode_type      = 'ka';
% A_basis_cmode   = 19.04;    %unitless frequency normalized by c/L

%Select mode index:
%This ensures that the gain center lies nearly on top of a cavity mode
cmode_type      = 'Na'; 
A_basis_cmode   = 20;         %must be an integer values

%Select refractive index: FP cavity only supports uniform index
A_basis_n       = 3 + 1i*0.0001;

%Select number of basis vectors, generated symmetrically around cmode frequency
A_basis_nCF     = 5;

%Specify \gamma_\perp:
A_cftd_gper     = 4;          %unitless frequency normalized by c/L

%Specify \gama_\parallel:
A_cftd_gpar     = 0.001;      %unitless frequency normalized by c/L

%Specify value of initial fields (usually a small number like 0.01)
A_cftd_eps      = 0.01;

% Add each array into the following cell of parameters. These parameters
% will be organized into param_vecs (below) according to the 'order'
% specified in the next line:
C_parameters = {A_basis_cmode.',A_cftd_gper.',A_cftd_gpar.',A_basis_n.',A_basis_nCF.',A_cftd_eps.'};

% Instruction for setup_genParameterVectors. Specifies order of output
% parameter vectors. Each number in this vector corresponds to a column
% number in param_vecs. For example, as configured below, values of
% \gamma_\perp will be stored in the fourth column of param_vecs and values
% of \gamma_\parallel will be stored in the fifth column. Instructions in
% setup_genParameterStruct expect this order.
order        = [1 4 5 2 3 6];

%Get param_vecs. Each row in param_vecs corresponds to parameters for
%exaclty one unique cavity based on the parameter value arrays
%specified above
param_vecs = setup_genParameterVectors(C_parameters,order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SET PUMP PARAMETERS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_pumpdata = struct('p_dens',10, ...    %density of pump per unit threshold
    't0',       0, ...                  %start time (unitless time normalized by L/c)
    't1',       40000, ...              %end time (unitless time normalized by L/c)
    'cratio',   0.1, ...                %calculation ratio (see description below)
    'sratio',   1, ...                  %save ratio (see description below)
    'prel_i',   1, ...                  %pump start value: as multiple of threshold
    'prel_f',   1.5, ...                %pump end value: as multiple of threshold
    'pg_len',   5, ...                  %length of each pump group (see description below)
    'xdens',    100, ...                %spatial points per period ([0,2*pi] range)
    'type',     'simple');              %pump type: 'simple' or 'hysteresis' (see description below)

% cratio (calculation ratio): t0 and t1 specify a time window of length
% T=t1-t0. cratio divides up the Runge-Kutta calculation into N segments of
% length T*cratio where N = ceil(1/cratio). For large calculations (ie:
% long T and large nCF), it is advisable to keep cratio small. Make sure
% that cratio*T is larger than at least a few round trips.

% sratio (save ratio): T*sratio specifies the fraction of the calculation that
% should be saved. T*sratio is a time duration at the very end of the time
% vector T and only values calculated in that duration are saved. We
% usually only need to analyze the calculation once it has reach steady
% state or periodic behavior, so saving only the last 'sratio' fraction of
% the calculation saves a lot of disk space.

% pg_len (pump group length): When jobs are being sent to the cluster, they
% are divided into J=K*G jobs where K is the number of rows in param_vecs 
% and G is determined by pg_len. For each kth job from the set of K jobs, 
% each job contains P=ceil((prel_f-prel_i)*p_dens) pump steps which are
% divided into groups of length pg_len. pg_len must contain at least one 
% pump step (ie. pg_len>=1). So for p_dens=10, and prel_i=1 and prel_f=1.5 
% we have P=5 total pump steps. If we choose pg_len=2, we need G= 3 
% pumpgroups to contain all pump steps. Then we have a total number of jobs 
% J=K*G=3 (the above parameters result in K=1).

% type ('simple' or 'hysteresis'): this is for pump "type". 'Simple'
% specifies that each pump step is independent. 'Hysteresis' seeds the next
% pump step with the values obtained from the end of the previous pump step
% calculation. In addition, 'hysteresis' starts at the lowest pump power
% and sequentially (or serially) goes to the highest pump power and then it
% goes backwards, seeding each lower pump step with values from the end
% of the calculation for the adjacent higher pump step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This generates the basis vectors and overlap integrals and stores all
%calculation details in S_setupdata and pump details in S_pumpdata.
%S_setupdata is a struct array of length K (same as param_vecs) where K is
%the number of unique cavities we will be simulating.
[S_setupdata,S_pumpdata] = setup_genParameterStruct(param_vecs,S_pumpdata,basis_type,cmode_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SET ODE SOLVER PARAMETERS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_odeoptions.benchmarking = true;       %to benchmark or not
S_odeoptions.use_ode_cpp  = true;       %true: use C++ Boost library. false: use MATLAB's ode45

% Only define the following if use_ode_cpp = true
S_odeoptions.log_level    = 'INFO';     %ERROR WARNING INFO DEGUB DEBUG1 DEBUG2 DEBUG3 DEBUG4 
S_odeoptions.odeint_const = true;       %true: use constant time step dt
S_odeoptions.dt           = 0.1;        %C++ constant time step duration
S_odeoptions.abs_error    = 1.0e-6;     %C++ ode absolute error
S_odeoptions.rel_error    = 1.0e-3;     %C++ ode relative error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Everything needed for calculations is stored in the following file
save('TD_init_parameters.mat','S_setupdata','S_pumpdata','param_vecs','S_odeoptions');

% Save instruction for number of pump groups as needed by bash 
% (only relevant when sending jobs to cluster)
f = fopen('for_bash','w');
fprintf(f,'%g',S_pumpdata.pgroups);
fclose(f);

%Create directory to store output from SLURM
mkdir('slurms');

%Output number of calulations and pump groups
leN = size(param_vecs,1);
fprintf(['Number of calculations: %g\n' ...
         'Number of pgroups: %g\n' ...
         'Length of pgroups: %g\n'], ...
         leN*S_pumpdata.pgroups,S_pumpdata.pgroups,S_pumpdata.pg_len);