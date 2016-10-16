function P_dot_x = get_P_dot_x(P_x_p,P_x_m,delta)
    P_dot_x   = (P_x_p - P_x_m)/(2*delta);
end