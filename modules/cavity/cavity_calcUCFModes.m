function [CFvals,CFvecs] = cavity_calcUCFModes(k_a,k_ref,nVec,nLeft,nRight,N_CF,nx,dx)
%possible values of type: openboth, openright, openleft

nsqrdVec = nVec.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(nx,1);
eigenmatrix = spdiags([-e 2*e -e], -1:1, nx, nx);
if nLeft~=inf
    eigenmatrix(1,1) = 2 - (2+1i*k_a*nLeft*dx)/(2-1i*k_a*nLeft*dx);
end
if nRight~=inf
    eigenmatrix(nx,nx) = 2 - (2+1i*k_a*nRight*dx)/(2-1i*k_a*nRight*dx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.tol=1e-9;
options.disp = 0;
% options.issym = 1;
inVindexmat = sparse(1:nx,1:nx,1./nsqrdVec);

if length(k_ref)==1
    tic;
    fprintf('Begin to solve the eigensystem.\n');
    [CFvecs, CFvals] = eigs(inVindexmat*eigenmatrix,N_CF,(k_ref*dx)^2,options);
    fprintf('Done.\n');
    toc;
    CFvals = sqrt(diag(CFvals))/dx;
elseif length(k_ref)==N_CF
    CFvals = zeros(N_CF,1);
    CFvecs = zeros(nx,N_CF);
    for ii=1:N_CF
        [CFvecs(:,ii), CFvals(ii)] = eigs(inVindexmat*eigenmatrix,1, (k_ref(ii)*dx)^2,options);
    end
    CFvals = sqrt(CFvals)/dx;
else
    fprintf('Error! The length of kref must be 1 or length(G_COM_M)\n');
    return;
end

for ii = 1:(N_CF)
    eta = sum( CFvecs(:,ii).^2.*nsqrdVec.')*dx;
    CFvecs(:,ii) = CFvecs(:,ii)/sqrt(eta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%