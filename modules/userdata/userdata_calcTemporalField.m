function field_t = userdata_calcTemporalField(coeffs,CFvecs,x,x0)
%     K = randperm(size(CFvecs,2))
%     [~,x_ind]   = max(abs(CFvecs(:,K(1))));
%     x_ind
    [~,x_ind] = min(abs(x - x0*x(end)));
    field_t   = squeeze(sum(coeffs.*repmat(CFvecs(x_ind,:),[size(coeffs,1),1,size(coeffs,3)]),2));
end