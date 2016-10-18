function dy = TDSSolvers_FP(t,y,g_per,g_par,ka,kt,A,D0_vec,D0,ncav)
    N = (length(A))^0.25;

    %Reshape inputs into identifiable vectors/arrays
    a_vec  = y(1:N);
    b_vec  = y(N+1:2*N);
    D_vec  = y(2*N+1:end);
    D_mat  = reshape(y((2*N+1):end),[N,N]);
    A_mat  = reshape(A,[N,N^3]);
    Ac_mat = A_mat;
    
    D0_vec = D0*D0_vec;
    NL_p = transpose(transpose(a_vec)*reshape(transpose(conj(b_vec))*A_mat,[N,N^2]));
    NL_m = transpose(transpose(conj(a_vec))*reshape(transpose(b_vec)*Ac_mat,[N,N^2]));

    dy = [(ka-kt.^2/ka)*1i*0.5.*a_vec + 1i*ka*0.5*b_vec/ncav^2; ...
          -g_per*b_vec - 1i*(a_vec.'*D_mat).'; ...
          -g_par*(D_vec - D0_vec) + 2*1i*(NL_p - NL_m)];
end

%Convert inputs to global
