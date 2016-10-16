function Dmn = get_Dmn_from_Dx(num,D_x,m,n)
    load('TD_init_parameters.mat','S_dparameters');
    dx = S_dparameters(num).dx;
    
    [~,CFvecs] = get_CF_basis(num);
    Dmn = sum(CFvecs(:,m).*D_x.*CFvecs(:,n))*dx;
end