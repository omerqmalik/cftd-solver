function Amat = salt_getAmat(CFvecs,dx,holeburning)
    Amat = transpose(CFvecs)*(CFvecs./repmat(1+holeburning,1,size(CFvecs,2)))*dx;
end