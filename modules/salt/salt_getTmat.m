function Tmat = salt_getTmat(lase_k,Amat,S_coredata,S_saltparameters)
    if strcmp(S_saltparameters.type,'SALT')
        Tmat = repmat(1i*lase_k^2*S_coredata.g_per/(S_coredata.g_per-1i*(lase_k-S_coredata.k_a))./(lase_k^2 - S_coredata.CFvals.^2),1,size(Amat,2)).*Amat;   %SALT
    elseif strcmp(S_saltparameters.type,'SVEA-SALT')
%         Tmat = repmat(1i*(S_coredata.k_a)^2*S_coredata.g_per/(S_coredata.g_per-1i*(lase_k-S_coredata.k_a))./((S_coredata.k_a^2 + 2*(lase_k-S_coredata.k_a)*S_coredata.k_a)-S_coredata.CFvals.^2),1,size(Amat,2)).*Amat;     %SVEA-SALT
%         Tmat = repmat(1i*(S_coredata.k_a)^2*S_coredata.g_per/(S_coredata.g_per-1i*lase_k)./(S_coredata.k_a^2 + 2*lase_k*S_coredata.k_a-S_coredata.CFvals.^2),1,size(Amat,2)).*Amat;     %SVEA-SALT
        Tmat = repmat(1i*(S_coredata.k_a)^2*S_coredata.g_per/(S_coredata.g_per+1i*(S_coredata.k_a-lase_k))./(2*lase_k*S_coredata.k_a - S_coredata.k_a^2 - S_coredata.CFvals.^2),1,size(Amat,2)).*Amat;
    end
end