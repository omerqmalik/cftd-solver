function noise_vec = get_noise_vec(N_CF,eps)
    noise_vec = zeros(1,N_CF^2+2*N_CF);
    noise_vec(1:2*N_CF) = eps;
    noise_vec((2:N_CF+1)*N_CF+(1:N_CF)) = eps;
end