function [Y_p,T] = get_Y_p_micro(num,D0_act)
    dirname  = ['micro/' num2str(num)];
    filename = ['p_' num2str(D0_act) '.mat'];
    load([dirname '/' filename]);
end