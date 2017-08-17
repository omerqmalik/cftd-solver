clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SET VECTOR PARAMETERS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basis_type = 'FP';

%for generating basis vectors
% A_basis_k_a   = 1.903995545007594e+03;        %unitless (k_a*L = w_a*L/c)
cmode_type = 'Na';
A_basis_cmode = 2000;        %unitless (k_a*L = w_a*L/c)
A_basis_n     = 3.3 + 0.0015*1i;
A_basis_nCF   = 7;

%for CFTD calculations
A_cftd_gper = 2.5;             %unitless (gper*L/c)
A_cftd_gpar = 1.25;    %unitless (gpar*L/c)
A_cftd_eps  = 0.01;

C_parameters = {A_basis_cmode.',A_cftd_gper.',A_cftd_gpar.',A_basis_n.',A_basis_nCF.',A_cftd_eps.'};
order        = [1 4 5 2 3 6];

param_vecs = setup_genParameterVectors(C_parameters,order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SET PUMP PARAMETERS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_pumpdata = struct('p_dens',5, ...   %density of pump per unit threshold
    't0',       0, ...                %start time (unitless: t*c/L)
    't1',       1000, ...             %end time (unitless: t*c/L)
    'cratio',   0.1, ...             %calculation ratio
    'sratio',   0.25, ...             %save ratio
    'prel_i',   1, ...               %as ratio of th
    'prel_f',   2, ...               %as ratio of th
    'pg_len',   20, ...               %length of each pump group
    'xdens',    100, ...                 %points per [0,2*pi] range
    'type',     'simple');          %pump type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S_setupdata,S_pumpdata] = setup_genParameterStruct(param_vecs,S_pumpdata,basis_type,cmode_type);
save('TD_init_parameters.mat','S_setupdata','S_pumpdata','param_vecs');

f = fopen('for_bash','w');
fprintf(f,'%g',S_pumpdata.pgroups);
fclose(f);

mkdir('slurms');

leN = size(param_vecs,1);
fprintf(['Number of calculations: %g\n' ...
         'Number of pgroups: %g\n' ...
         'Length of pgroups: %g\n'], ...
         leN*S_pumpdata.pgroups,S_pumpdata.pgroups,S_pumpdata.pg_len);
