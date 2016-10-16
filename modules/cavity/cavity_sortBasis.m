function [CFvals,CFvecs,m] = cavity_sortBasis(CFvals,CFvecs,k_a)
    midcf = ceil(size(CFvecs,3)/2);
    CFvecs = CFvecs(:,:,midcf);
    CFvals = CFvals(:,midcf);
    [~,m] = sort(abs(real(CFvals) - k_a));
    CFvals = CFvals(m);
    CFvecs = CFvecs(:,m);
end
