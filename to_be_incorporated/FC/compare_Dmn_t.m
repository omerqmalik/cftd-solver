function compare_Dmn_t(num,D0_act,t0,t1,skip,m,n)
    figure;
    load('TD_init_parameters.mat','S_dparameters');
    nCF = S_dparameters(num).nCF;
    
    [Dmn,t] = get_Dmn_micro(num,D0_act);
    [~,ind0] = get_closest_match(t0,t);
    [~,ind1] = get_closest_match(t1,t);
    
    inds = ind0:skip:ind1;
    t_this = t(inds);
    mn_ind = nCF*(m-1) + n;
    Dmn_orig = Dmn(inds,mn_ind);
    
    plot(t_this,Dmn_orig);
    hold on;
    
    t_len = length(t_this);
    Dmn_calc = zeros(1,t_len);
    for i = 1:length(t_this)
        D_x_this = get_D_x_sg(num,D0_act,t_this(i));
        Dmn_calc(i) = get_Dmn_from_Dx(num,D_x_this,m,n);
    end
    
    plot(t_this,real(Dmn_calc));
end