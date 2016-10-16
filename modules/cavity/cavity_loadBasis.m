function [CFvals,CFvecs,basis_type] = cavity_loadBasis(basis_loc,k_a)
    load(basis_loc,'CFvals','CFvecs','basis_type');
    [CFvals,CFvecs] = cavity_sortBasis(CFvals,CFvecs,k_a);
end
