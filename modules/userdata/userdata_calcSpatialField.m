function field_x = userdata_calcSpatialField(coeffs,CFvecs)
    field_x = sum(CFvecs.*repmat(coeffs,[size(CFvecs,1),1]),2);
end