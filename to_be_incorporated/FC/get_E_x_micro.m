function [E_x_CF,x,dt] = get_E_x_micro(num,D0_act,t0,inds)
	load('TD_init_parameters.mat','S_dparameters');
    x = S_dparameters(num).x;
    nCF = S_dparameters(num).nCF;

    dirname  = ['micro/' num2str(num)];
    filename = ['e_' num2str(D0_act) '.mat'];
    load([dirname '/' filename]);
    
    t = T;
    Y = Y_e;
    
    [~,this_CFvecs] = get_CF_basis(num);
    this_CFvecs = this_CFvecs(:,1:size(Y,2));
    
    ensL = length(inds);
    [~,t_ind] = get_closest_match(t0,t);
    
    E_x = zeros(size(this_CFvecs,1),ensL);
    for i = 1:ensL
        ind_this = t_ind + inds(i);
        E_x_CF = zeros(size(this_CFvecs));
        for ii = 1:nCF
            E_x_CF(:,ii) = Y(ind_this,ii)*this_CFvecs(:,ii);
        end
        E_x_CF = sum(E_x_CF,2);
        E_x(:,i) = E_x_CF;
    end
    dt = t(2) - t(1);
end