function field_t = userdata_calcTemporalField(coeffs,CFvecs)    
    [~,x_ind]   = max(abs(CFvecs(:,1)));
    field_t = squeeze(sum(coeffs.*repmat(CFvecs(x_ind,:),[size(coeffs,1),1,size(coeffs,3)]),2));
end