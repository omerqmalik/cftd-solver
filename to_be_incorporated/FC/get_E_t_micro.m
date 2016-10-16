function [E_t,T] = get_E_t_micro(num,D0_act)
    dirname  = ['micro/' num2str(num)];
    filename = ['E_t_' num2str(D0_act) '.mat'];
    load([dirname '/' filename]);
end