function [D_x,x] = get_D_x_sg(num,D0_act,t0)
    %load basic parameters
    load TD_init_parameters;
    g_per     = S_dparameters(num).g_per;
    
    %get E_x
    [E_x,x,dt] = get_E_x_micro(num,D0_act,t0,0);
    
    %get P_x
    P_x = get_P_x_micro(num,D0_act,t0,[-1 0 1]);
    
    %Calculate D_t_sg
    D_x = get_D_x_sg_from_fields(E_x,P_x,dt,g_per);
end

function D_x = get_D_x_sg_from_fields(E_x,P_x,dt,g_per)
    %Get P_dot and then calculate D_t + SG smoothing
    P_dot_x = get_P_dot_x(P_x(:,3),P_x(:,1),dt);
    D_x = get_D_t(P_dot_x,P_x(:,2),E_x,g_per);
end