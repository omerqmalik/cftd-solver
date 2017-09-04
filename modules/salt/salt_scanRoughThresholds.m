function S_rough = salt_scanRoughThresholds(S_coredata,S_saltparameters) %(CcfvasFInt,CcfvesFInt,fn_tuning)

TN_skCtr = S_saltparameters.skCtr;
TN_skWnd = S_saltparameters.skWnd;
TN_skN   = S_saltparameters.skN;
TN_sskN  = S_saltparameters.sskN;

s_kd  = TN_skWnd/TN_skN;          %Length (d) of each slice of k-values
ss_kd = s_kd/TN_sskN;             %Length (d) of each subslice of k-values

[~,CFvecs,~,dx] = setup_getBasis(S_coredata.basis_loc);
Amat = salt_getAmat(CFvecs,dx,zeros(size(CFvecs,1),1));      %Needed to calculate threshold matrix

rghsol_N = 0;               %No rough solutions have been found yet
for s_i = 1:TN_skN             %For each slice
    fprintf('slice # %g of %g\n',s_i,TN_skN);
    
    knum_i = (TN_skCtr-TN_skWnd/2)+(s_i-0.5)*s_kd;
    
    for ss_i = 1:TN_sskN                %For each subslice
        rghlase_k = knum_i - s_kd/2 + ss_kd*ss_i;
        Tmat = salt_getTmat(rghlase_k,Amat,S_coredata,S_saltparameters);       %Calculate threshold matrix
        [tves_this,tvas_this] = eig(Tmat);                  %Solve threshold matrix eigenvalue problem
        tvas_this = nonzeros(tvas_this);                    %Filter out the nonzero elements
        
        if (ss_i==1)
            tvas_prev = tvas_this;
        else
            [is_sol,thSol,tvesSol,tvas_prev] = hlpr_findCrossoverSol(tvas_this,tves_this,tvas_prev);
            if is_sol==1        %Crossover found for this particular subslice
                rghsol_N = rghsol_N + 1;
                kRough(rghsol_N) = rghlase_k;
                thRough(rghsol_N) = thSol;
                tvesRough(:,rghsol_N) = tvesSol;
                fprintf('\nRough Solution #%g found!\nsubslice = %g\nkRough = %g\nthRough = %g\n\n',rghsol_N,ss_i,kRough(rghsol_N),thRough(rghsol_N));
            end
        end
    end
end
%Keep only positive values
kRough              = kRough(thRough>0);
tvesRough           = tvesRough(:,thRough>0);
thRough             = thRough(thRough>0);

%Sort
[S_rough.th,S_rough.ord] = sort(thRough,'ascend');
S_rough.k                = kRough(S_rough.ord);
S_rough.tves             = tvesRough(:,S_rough.ord);

fprintf('\n');
end