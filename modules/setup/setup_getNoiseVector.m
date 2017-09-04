function noise_vec = setup_getNoiseVector(nCF,eps)
    noise_vec = zeros(nCF^2+2*nCF,1);
    noise_vec(1:2*nCF) = eps;
    noise_vec((2:nCF+1)*nCF+(1:nCF)) = eps;
end
