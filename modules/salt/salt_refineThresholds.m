function S_fine = salt_refineThresholds(S_rough,S_coredata,S_saltparameters,holeburning)

TN_ntest  = S_saltparameters.ntest;
TN_therr  = S_saltparameters.therr;
TN_kerr   = S_saltparameters.kerr;

kRough    = S_rough.k;
thRough   = S_rough.th;
tvesRough = S_rough.tves;
ordRough  = S_rough.ord;

[~,CFvecs,~,dx] = setup_getBasis(S_coredata.basis_loc);
Amat = salt_getAmat(CFvecs,dx,holeburning);

% Improve the accuracy of the first few thresholds

for m = 1:length(kRough) % Calculate the threshold for the "mth" mode
    fprintf('Refining rough mode #: %g\n',m);
    rghlase_k = kRough(m);
    rghlase_th = thRough(m);

    for i = 1:TN_ntest          %keep iterating until if clause validates convergence
        fprintf('iteration #: %g\n',i);
        [knumTmp,thTmp,tvesFine(:,m)] = salt_scanFineThresholds(rghlase_k,Amat,tvesRough(:,m),S_coredata,S_saltparameters);
        % added by hand to catch the slow-convergence one
        fprintf('%10.10f,%10.10f\n',abs(rghlase_th-thTmp),abs(rghlase_k-knumTmp));
        if (abs(rghlase_th-thTmp)<TN_therr && abs(rghlase_k-knumTmp)<TN_kerr)
            kFine(m) = knumTmp;
            thFine(m) = rghlase_th;
            maxpole_pos(m) = find(abs(tvesFine(:,m))==max(abs(tvesFine(:,m))));     %Find the maximum amplitude pole
            
            fprintf('Fine solution is good enough. Done!\n\n\n');
            break;  %convergence validated, solution is fine enough
        else
            fprintf('Fine solution not good enough. Keep iterating!\n\n');
        end
        rghlase_k = knumTmp;
        rghlase_th = thTmp;
    end
end

%Sorting and ordering
[S_fine.th,ord]   = sort(thFine);
S_fine.k          = kFine(ord);
S_fine.tves       = tvesFine(:,ord);
S_fine.maxpolepos = maxpole_pos(ord);
S_fine.ord        = ordRough(ord);

for i=length(kRough):-1:2
    if norm(S_fine.tves(:,i)-S_fine.tves(:,i-1))<1e-4  && abs(S_fine.k(i)-S_fine.k(i-1))<1e-4       %If true, then consecutive quantities belong to the same solution
        S_fine.k(i)=[];
        S_fine.th(i)=[];
        S_fine.tves(:,i)=[];
        S_fine.maxpolepos(:,i)=[];
    end
end

fprintf('Threshold solver successful!\n\n\n');