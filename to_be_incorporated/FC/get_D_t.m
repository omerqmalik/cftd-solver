function D_t = get_D_t(P_dot,P_t,E_t,g_per)
    D_t = 1i*(P_dot + g_per*P_t)./(g_per*E_t);
end