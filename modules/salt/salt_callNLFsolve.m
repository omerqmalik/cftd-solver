function [S_nlsol,holeburning,exitflag] = salt_callNLFsolve(d0,on_modes,S_crit,S_coredata,S_saltparameters)
%M: total number of modes
%Todo: make functions of intermediate parts

for i=1:on_modes
    amax_this                        = S_crit.a(S_crit.maxpolepos(i),i);   %the max element of the current "direction" vector
    S_crit.a(:,i)                    = S_crit.a(:,i)/amax_this;            %This prevents convergence to the trivial solution (see Li's thesis, Ch. 4, pg. 78)
    S_crit.a(S_crit.maxpolepos(i),i) = amax_this + 1i*S_crit.k(i);     
end
a_rshpd = reshape(S_crit.a,S_coredata.nCF*on_modes,1);     %a contains direction vectors (which consist of cf-state coefficients for each of the G_COM_M cf-states). Reshape a so that they're all in a single column
a_RI = [imag(a_rshpd);real(a_rshpd)];    %a vector separated into Real and Imaginary components

%Find the "direction" vector for each mode that has turned ON, acocunt for nonlinear interactions
[~,~,~,dx] = setup_getBasis(S_coredata.basis_loc);
options=optimset('Display','iter','TolFun',S_saltparameters.MMfslvOptTolFun,'TolX',S_saltparameters.MMfslvOptTolX,'MaxIter',500,'MaxFunEvals',4000,'PlotFcn',@optimplotfirstorderopt,'Algorithm','levenberg-marquardt');
[agsol,~,exitflag,~] = fsolve(@(input)salt_calcSALTEquation(input,d0,on_modes,S_crit,S_coredata,S_saltparameters,dx),a_RI,options);

if exitflag < 1         %Check here whether the nonlinear solver was able to find a non-trivial solution or not
    S_nlsol.a = [];
    S_nlsol.k = [];
    holeburning = [];
    
    fprintf('Uh-oh, fsolve did not converge. Exitflag: %d\n', exitflag);
    
    return;
end

agsol = reshape(agsol,S_coredata.nCF,2*on_modes);
S_nlsol.k = zeros(1,on_modes);
S_nlsol.a = zeros(S_coredata.nCF,on_modes);
for i=on_modes:-1:1
    [S_nlsol.k(i),S_nlsol.a(:,i)] = get_Normalized_asol(agsol,S_crit.maxpolepos,on_modes,i,S_coredata.nCF);        %No longer normalized nor gauge solution
    
    if max(abs(S_nlsol.a(:,i))) < 1e-12       %If the solution is too small
        fprintf('Warning: Trivial Solution!');
        return;
    end
end

holeburning = calc_holeburning(on_modes,S_nlsol.k,S_nlsol.a,S_coredata);
end

function holeburning = calc_holeburning(on_modes,ksol,asol,S_coredata)
    holeburning      = zeros(size(S_coredata.CFvecs,1),1);

    for i=1:on_modes
        holeburning = holeburning + abs(S_coredata.CFvecs*asol(:,i)).^2/(1+((ksol(i)-S_coredata.k_a)/S_coredata.g_per)^2);
    end
end

function [ksol,asol] = get_Normalized_asol(sol,maxpos,M,mode,nCF)
    ksol = sol(maxpos(mode),mode);      %The lasing k solution for this mode
    sol(maxpos(mode),mode)=0;           %Set imaginary part of fsolve solution to zero at max pole for this mode
    asol(1:nCF) = 1i*sol(:,mode)+sol(:,mode+M);   %Reconstruct the a-vector for this mode
    amax = asol(maxpos(mode));        %Get value at max pole of a-vector for this mode
    asol = asol*amax;                %Why do we do this? To undo the effect of normalization. Go from bmu to amu
    asol(maxpos(mode)) = amax;        %Normalize
end