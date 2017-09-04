function [kOut,thOut,tvesOut] = salt_scanFineThresholds(knum,Amat,tvesRough,S_coredata,S_saltparameters)

TN_fskWnd = S_saltparameters.fskWnd;
TN_fskN   = S_saltparameters.fskN;
TN_kcri   = S_saltparameters.kcri;
TN_kprec  = S_saltparameters.kprec;
TN_dkprec = S_saltparameters.dkprec;

N_totalCF = S_coredata.nCF;

kFiner     = zeros(1,N_totalCF);      %A maximum of G_COM_M eigenvalues could have crossed over
thFiner    = zeros(1,N_totalCF);
tvesFiner  = zeros(N_totalCF,N_totalCF);
finersol_N = 0;                 %No fine solutions found just yet

fs_kCtr = knum+TN_fskWnd/2;       %Center of fine-slice window
fs_kd = TN_fskWnd/TN_fskN;          %fine-slice step size

while(1)
    for fs_i = 1:TN_fskN  %For each fineslice
        %fprintf('subiteration %g of %g\n',fs_i,fs_kN)
        finelase_k = (fs_kCtr+TN_fskWnd/2)-fs_i*fs_kd;           %This time go backwards, need to crossover in the other direction
        Tmat = salt_getTmat(finelase_k,Amat,S_coredata,S_saltparameters);         %Calculate threshold matrix
        [tves_this,tvas_this] = eig(Tmat);                     %This segment is identical to finding the rough thresholds
        tvas_this = nonzeros(tvas_this);
        if (fs_i==1)
            tvas_prev = tvas_this;
        else
            [is_sol,thSol,tvesSol,tvas_prev] = salt_findCrossoverSol(tvas_this,tves_this,tvas_prev);
            if is_sol==1        %If a fine solution is found
                finersol_N = finersol_N + 1;        %increase solution count
                kFiner(finersol_N) = finelase_k;    %record finer k value
                thFiner(finersol_N) = thSol;        %record finer threshold
                tvesFiner(:,finersol_N) = tvesSol;  %record finer 'direction'
            end
        end
    end
    
    if(finersol_N==0)               %If no finer solution is found (yet)
        TN_fskWnd = TN_fskWnd*2;    %Increase search window
        fs_kd = fs_kd*2;            %Larger window needs more slices
        continue;                   %Try again with larger window
    else                            %Finer solutions found!
        [err,p] = howCloseToRoughSoln(finersol_N,tvesFiner,tvesRough);
        if err > 2                  %If the solution isn't fine enough (Why 2?)
            finersol_N=0;           %Try again with larger window
            TN_fskWnd = TN_fskWnd*2;
            fs_kd = fs_kd*2;
            continue;               %Try again
        end
        
        kTmp = kFiner(p);               %Might need fine tuning
        thTmp = thFiner(p);             %Might need fine tuning
        tvesOut = tvesFiner(:,p);       %Good enough, let's keep it (Why?)
        if ( abs(kTmp-knum) > TN_kcri)     %Change is still too large
            kOut = kTmp;                %Use this k value for now
            thOut = thTmp;              %Use this threshold for now
            return;             %Now go back and start with kOut/thOut and repeat process
        end

        if ( abs(kTmp-knum) < TN_kprec || fs_kd < TN_dkprec )     %It's fine enough (solution found!)
            kOut = kTmp;    %Final solution
            thOut = thTmp;  %Final threshold
            fprintf('Fine solution found!\nlasingk = %6.4f\nD0 = %6.4f \n',knum,thOut);
            return;     %We're done here
        else
            %Get ready for another iteration...
            knum = kTmp;                %Move the search window to kTmp
            TN_fskWnd =TN_fskWnd/10;        %Try a smaller window
            fs_kd = 2*TN_fskWnd/TN_fskN;
            finersol_N = 0;
            kFiner = zeros(1,N_totalCF);
            thFiner = zeros(1,N_totalCF);
            tvesFiner = zeros(N_totalCF,N_totalCF);
        end
    end

end
end

function [err_min,p] = howCloseToRoughSoln(finersol_N,tvesFiner,tvesRough)
errorsum = zeros(1,length(finersol_N));
for i=1:finersol_N
    errorsum(i) = sum(abs(abs(tvesFiner(:,i))-abs(tvesRough)));     %Best solution should be close to the rough solution. This ascertains that the solution hasn't moved too far in the *other* direction
end
[err_min,p] = min(errorsum);        %Find the closest one
end