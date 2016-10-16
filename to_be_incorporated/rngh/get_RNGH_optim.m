function RNGH_optim = get_RNGH_optim(kappa,alpha,lambda_guess,gpar_guess)
    options = optimoptions('fsolve','Display','final-detailed','MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-18,'TolX',1e-18);
    
    eflag = 0;
    while eflag < 1
        [RNGH_u,fval,eflag] = fsolve(@(x) RNGH_root_u(x,kappa,alpha),[lambda_guess gpar_guess],options);
        gpar_guess = gpar_guess/2;
    end
    
    RNGH_optim = RNGH_u;
end

function y = RNGH_root_u(x,kappa,alpha)
    lambda = x(1);
    gpar   = x(2);
    
    eqn  = get_RNGH_upper_bound(lambda,gpar,kappa);
    y(1) = alpha - real(eqn);    y(2) = imag(eqn);
end