function [Dmn,t] = get_Dmn_micro(num,D0_act)
	load('TD_init_parameters.mat','S_dparameters');

    dirname  = ['micro/' num2str(num)];
    filename = ['d_' num2str(D0_act) '.mat'];
    load([dirname '/' filename]);
    
    t = T;
    Dmn = Y_d;
end