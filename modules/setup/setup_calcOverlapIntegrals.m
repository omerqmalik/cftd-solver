function [A,A_index,B,D0_vec,D_index] = setup_calcOverlapIntegrals(nCF,n_sqrd,pump_profile,CFvecs,basis_type,dx,M)
    A  = zeros(1,nCF^4);
    B  = zeros(1,nCF^2);
    
    D0_vec  = zeros(nCF^2,1);
    A_index = zeros(nCF^4,4);
    D_index = zeros(nCF^2,2);
  
    tic;
    fprintf('Begin calculation of overlap integrals.\n');
    for i = 1:nCF     %m
        fprintf('1: %g\n',i);
        for j = 1:nCF     %n
            for ii = 1:nCF    %r
                for jj = 1:nCF    %s
                    if strcmp(basis_type,'UCF')
                        A(jj + (ii-1)*nCF + nCF^2*(j-1) + nCF^3*(i-1)) = sum(n_sqrd.*CFvecs(:,i).*CFvecs(:,ii).*conj(CFvecs(:,jj)).*CFvecs(:,j))*dx; %A_mrsn
                    elseif strcmp(basis_type,'FP')
                        m = M(i); r = M(ii); s = M(jj); n = M(j);
                        A(jj + (ii-1)*nCF + nCF^2*(j-1) + nCF^3*(i-1)) = (-(m-r-s-n==0) + (m+r-s-n==0) + (m-r+s-n==0) - (m+r+s-n==0) + (m-r-s+n==0) - (m+r-s+n==0) - (m-r+s+n==0))/2;
                    elseif strcmp(basis_type,'RING')
                        m = M(i); r = M(ii); s = M(jj); n = M(j);
                        A(jj + (ii-1)*nCF + nCF^2*(j-1) + nCF^3*(i-1)) = (-m+r-s+n==0); %A_mrsn
                    end
                    A_index(jj + (ii-1)*nCF + nCF^2*(j-1) + nCF^3*(i-1),:) = [i,ii,jj,j];
                end
            end
        end
    end
    A = reshape(A,[nCF^2,nCF^2]).';
    A_index = permute(reshape(A_index,[nCF^2,nCF^2,4]),[2 1 3]);
    fprintf('Done.\n');
    toc;

    for i = 1:nCF     %n
        fprintf('2: %g\n',i);
        for j = 1:nCF     %m
            if strcmp(basis_type,'UCF')
                B(j + (i-1)*nCF) = sum(CFvecs(:,i).*CFvecs(:,j))*dx;
            end
        end
    end
    B = reshape(B,[nCF,nCF]);

    for i = 1:nCF
        fprintf('3: %g\n',i);
        for j = 1:nCF
            if strcmp(basis_type,'UCF')
                D0_vec(j + (i-1)*nCF) = sum(n_sqrd.*CFvecs(:,i).*pump_profile.*CFvecs(:,j))*dx;
            elseif strcmp(basis_type,'FP')
                D0_vec(j + (i-1)*nCF) = sum(CFvecs(:,i).*pump_profile.*CFvecs(:,j))*dx;
            elseif strcmp(basis_type,'RING')
                D0_vec(j + (i-1)*nCF) = sum(conj(CFvecs(:,i)).*pump_profile.*CFvecs(:,j))*dx;
            end
%             fprintf('%g%g\n',j,i);
            D_index(j + (i-1)*nCF,:) = [i,j];
        end
    end
    D_index = reshape(D_index,[nCF,nCF,2]);
end

%     An = reshape(permute(reshape(An,[nCF nCF nCF nCF]),[4 2 1 3]),[nCF^2,nCF^2]);
%     An_index = permute(reshape(An_index(:,[4 2 1 3]),[nCF^2,nCF^2,4]),[2 1 3]);