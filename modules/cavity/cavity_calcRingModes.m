function [CFvals,CFvecs,dx,nx,x,w_FSR,Na,M] = cavity_calcRingModes(k_a,n,nCF,xdens)
    %units: 
    %k_a    (unitless, k_a*L)
    %w_FSR  (unitless, w_FSR*L/c)
    %CFvals (unitless, CFvals*L/c)
    %CFvecs (unitless, by definition)
    %x,dx   (unitless, x/L)
    
    w_FSR = 2*pi/n;                     %FSR
    Na = round(k_a/real(w_FSR));        %central mode number
    if mod(nCF,2) == 0
        M = (Na-nCF/2):(Na+(nCF/2-1));
    else
        M  = (Na-(nCF-1)/2):(Na+(nCF-1)/2); %modes to include
    end
    
    dx = 2/max(M)/xdens;
    nx = 1/dx;
    x = linspace(0,dx*nx,nx);
    
    CFvals = zeros(nCF,1);
    CFvecs = zeros(length(x),nCF);
    index = 1;
    for m = M
        CFvals(index)  = w_FSR*m;
        CFvecs(:,index) = exp(1i*2*pi*m*x);
        index = index + 1;
    end
    
    [CFvals,CFvecs,m] = cavity_sortBasis(CFvals,CFvecs,k_a);
    M = M(m);
end
