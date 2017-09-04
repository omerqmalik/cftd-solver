function dy = MBSolvers_inversion(t,y,pump,g_par,Efield,Pfield,field_t)
    E = interp1(field_t,Efield,t);
    P = interp1(field_t,Pfield,t);
    
    %Note: The 'pump' is assumed to be uniform here
    dy = g_par*(pump - y) - g_par*imag(E*conj(P));
end