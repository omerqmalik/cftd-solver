function field_x = userdata_calcSpatialField(coeffs,CFvecs)
%     field_x = sum(CFvecs.*repmat(coeffs,[size(CFvecs,1),1]),2);
    field_x = reshape(sum(repmat(CFvecs,[size(coeffs,1),1]).*repelem(coeffs,size(CFvecs,1),1),2),[size(CFvecs,1),size(coeffs,1)]);
end