function dy = TDSSolvers_UCF(t,y,g_per,g_par,k_a,CFvals,A,B,pump_pwr,N)
    %Reshape inputs into identifiable vectors/arrays
    a_vec  = y(1:N);
    b_vec  = y(N+1:2*N);
    D_vec  = y((2*N+1):end);
    D_mat  = reshape(D_vec,[N,N]);

    %Consult ring laser solver
    Q = conj(b_vec)*a_vec.';
    NL_term = A*reshape(Q-Q',[N^2,1]);

    dy = [(k_a-CFvals.^2/k_a)*1i*0.5.*a_vec + 1i*k_a*0.5*(b_vec.'*B).'; ...
          -g_per*b_vec - 1i*g_per*(a_vec.'*D_mat).'; ...
          -g_par*(D_vec - pump_pwr) + 1i*g_par*NL_term];
end