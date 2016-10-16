function Dth = cavity_calcDth(n,gper,nu0,kappa0,k_a)
    delta0    = -gper*(k_a^2 - nu0^2 + kappa0^2)/2/(gper*k_a + nu0*kappa0);
    Dth       = real(n)^2*(2*gper*nu0*kappa0 - 2*k_a*delta0^2 - (k_a^2 - nu0^2 + kappa0^2)*delta0)/gper/k_a^2;
end
