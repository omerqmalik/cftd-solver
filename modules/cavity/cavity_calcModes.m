function [CFvals,CFvecs,dx,nx,x,w_FSR,Na,M] = cavity_calcModes(basis_type,k_a,n,nCF,xdens)
    %units: 
    %k_a    (unitless, k_a*L)
    %w_FSR  (unitless, w_FSR*L/c)
    %CFvals (unitless, CFvals*L/c)
    %CFvecs (unitless, by definition)
    %x,dx   (unitless, x/L)
    
    if strcmp(basis_type,'RING')
        w_FSR = 2*pi/n;
        basis_func = @(m,x) exp(1i*2*pi*m*x);
    elseif strcmp(basis_type,'FP') || strcmp(basis_type,'UCF')   %Use FP basis as an estimate of w_FSR and basis_func
        w_FSR = pi/n;
        basis_func = @(m,x) sqrt(2)*sin(pi*m*x);
    end
    
    Na = round(k_a/real(w_FSR));        %central mode number
    if mod(nCF,2) == 0
        M = (Na-nCF/2):(Na+(nCF/2-1));
    else
        M  = (Na-(nCF-1)/2):(Na+(nCF-1)/2); %modes to include
    end
    
    dx = 2/max(M)/xdens;
    nx = 1/dx;
    x = linspace(0,dx*nx,fix(nx));
    
    CFvals = zeros(nCF,1);
    CFvecs = zeros(length(x),nCF);
    
    if strcmp(basis_type,'RING') || strcmp(basis_type,'FP')
        index = 1;
        for m = M
            CFvals(index)  = w_FSR*m;
            CFvecs(:,index) = basis_func(m,x);
            index = index + 1;
        end
    elseif strcmp(basis_type,'UCF')
        nVec = ones(1,length(x))*n;
        [CFvals,CFvecs] = cavity_calcUCFModes(k_a,k_a,nVec,1,1,nCF,fix(nx),dx);
        w_FSR = mean(diff(sort(real(CFvals))));
    end
    
    [CFvals,CFvecs,m] = cavity_sortBasis(CFvals,CFvecs,k_a);
    M = M(m);
end