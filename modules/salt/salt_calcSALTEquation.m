function fval=salt_calcSALTEquation(input,d0,on_modes,S_crit,S_coredata,S_saltparameters,dx)
% FVAL=TS(INPUT,MODENUM) is the nonlinear TS equation
% Here use normalized a_i = a_i/a_dominant.
    
    holeburning = zeros(size(S_coredata.CFvecs,1),1);                %The nonlinear term (interaction component)

    a_gauge = zeros(S_coredata.nCF,on_modes);             %Holds "direction" vector with gauge specified
    amax_gauge = zeros(1,on_modes);                %Holds max of "direction" vector with gauge specified
    fval = zeros(S_coredata.nCF,1);
    lase_k = zeros(S_coredata.nCF,1);              %The fine-tuned output lasing frequency

    %Set gauge condition of "direction" vector and calculate non-linear term
    atmp = reshape(input,S_coredata.nCF,2*on_modes);
    for i=1:on_modes
        [a_gauge(:,i),amax_gauge(i),lase_k(i)] = get_agauge(atmp,S_crit.maxpolepos,on_modes,i,S_coredata.nCF);      %apply gauge condition
        %Calculation of non-linear term (Eq 4.9 Li's Thesis)
        holeburning = holeburning + amax_gauge(i)^2*abs(S_coredata.CFvecs*a_gauge(:,i)).^2/(1+((lase_k(i)-S_coredata.k_a)/S_coredata.g_per)^2);
    end

    %Calculate all terms of the non-linear SALT equation
    for i=1:on_modes               %For each mode below critical value (ie: current ON modes)
        for j=1:S_coredata.nCF     %For every eigenvalue
            modal_interaction = sum(S_coredata.CFvecs(:,j).*S_coredata.CFvecs*a_gauge(:,i)./(1 + holeburning))*dx;
            
            if strcmp(S_saltparameters.type,'SALT')
                tmp = 1i*d0*lase_k(i)^2*S_coredata.g_per/(S_coredata.g_per-1i*(lase_k(i)-S_coredata.k_a))/(lase_k(i)^2-S_coredata.CFvals(j)^2)*modal_interaction;   %SALT
            elseif strcmp(S_saltparameters.type,'SVEA-SALT')
%                 tmp = 1i*d0*(S_coredata.k_a^2)*S_coredata.g_per/(S_coredata.g_per-1i*(lase_k(i)-S_coredata.k_a))/((S_coredata.k_a^2 + 2*(lase_k(i)-S_coredata.k_a)*S_coredata.k_a)-S_coredata.CFvals(j)^2)*modal_interaction; %SVEA-SALT
                tmp = 1i*d0*(S_coredata.k_a^2)*S_coredata.g_per/(S_coredata.g_per+1i*(S_coredata.k_a-lase_k(i)))/(2*lase_k(i)*S_coredata.k_a-S_coredata.k_a^2-S_coredata.CFvals(j)^2)*modal_interaction; %SVEA-SALT
            end

            fval((i-1)*2*S_coredata.nCF+2*j-1) = real(tmp-a_gauge(j,i));   %Put real part in single vector
            fval((i-1)*2*S_coredata.nCF+2*j) = imag(tmp-a_gauge(j,i));     %Put imag part in (same) single vector
        end
    end
end

function [ag,agmax,k] = get_agauge(atmp,maxpos,M,mode,nCF)
    k = atmp(maxpos(mode),mode);           %Extract lasing k from the "direction" vector
    atmp(maxpos(mode),mode) = 0;           %Set imaginary component to zero (gauge condition)
    ag(1:nCF) = atmp(:,mode+M)+ 1i*atmp(:,mode);     %Reconstruct matrix containing a-vectors wth gauge condition applied
    agmax = ag(maxpos(mode));              %Record the maximum value of a-vector before normalization
    ag(maxpos(mode)) = 1;                  %Normalize maximum value
end