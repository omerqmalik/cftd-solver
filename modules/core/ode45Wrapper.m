function [T, Y] = ode45Wrapper(S_coredata, pstep, tv ,noise_vec, opts)
%ODE45WRAPPER Summary of this function goes here
%   Detailed explanation goes here

global usingc;

    g_per      = S_coredata.g_per;
    g_par      = S_coredata.g_par;
    k_a        = S_coredata.k_a;
    CFvals     = S_coredata.CFvals;
    A          = S_coredata.A;
    B          = S_coredata.B;
    pump       = S_coredata.pump;
    D0_vec     = S_coredata.D0_vec;
    n          = S_coredata.n;
    basis_type = S_coredata.basis_type;
    pump_pwr = pump(pstep)*D0_vec;

if (usingc)
    
    
    t_initial = tv(1);
    t_final = tv(2);
    
    [T, Y] = runge_kutta4(S_coredata, pump_pwr, t_initial, t_final, noise_vec, basis_type);
    
else
    
    
    if strcmp(basis_type,'RING')
        solver_func = @(t,y) TDSSolvers_RING(t,y,g_per,g_par,k_a,CFvals,A,pump_pwr,length(CFvals),n);
    elseif strcmp(basis_type,'FP')
        solver_func = @(t,y) TDSSolvers_FP(t,y,g_per,g_par,k_a,CFvals,A,pump_pwr,length(CFvals),n);
    elseif strcmp(basis_type,'UCF')
        solver_func = @(t,y) TDSSolvers_UCF(t,y,g_per,g_par,k_a,CFvals,A,B,pump_pwr,length(CFvals));
    end

    [T, Y] = ode45(solver_func, tv ,noise_vec, opts);

end


end

