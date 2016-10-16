function [T_out,P_dot,P_t_out,E_t_out] = get_P_dot_t(T_in,P_t_in,E_t_in)
    P_t_p   = P_t_in(3:end);
    P_t_m   = P_t_in(1:(end-2));
    
    P_t_out = P_t_in(2:(end-1));
    E_t_out = E_t_in(2:(end-1));
    T_out   = T_in(2:(end-1));
    
    T_diff  = T_in(2) - T_in(1);
    P_dot   = (P_t_p - P_t_m)/(2*T_diff);
end