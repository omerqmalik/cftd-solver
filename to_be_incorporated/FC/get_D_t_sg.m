function [D_t_sg,D_t,T] = get_D_t_sg(num,D0_act,cratio,isplot,issmooth)
    %load basic parameters
    load TD_init_parameters;
    basis_loc = S_dparameters(num).basis_loc;
    k_a       = S_dparameters(num).k_a;
    g_per     = S_dparameters(num).g_per;
    
    %load E_t
    [E_t,T] = get_E_t_micro(num,D0_act);
    
    %calculate P_t
    Y_p = get_Y_p_micro(num,D0_act);
    P_t = get_E_t(Y_p,num,basis_loc,k_a);
    
    %Truncate E_t, P_t, T
    E_t = E_t(ceil((1-cratio)*length(E_t))+1:end);
    P_t = P_t(ceil((1-cratio)*length(P_t))+1:end);
    T   = T(ceil((1-cratio)*length(T))+1:end);
    
    %Calculate D_t_sg
    [D_t_sg,D_t,T] = get_D_t_sg_from_fields(T,E_t,P_t,g_per,isplot,issmooth);
end

function [D_t_sg,D_t,T] = get_D_t_sg_from_fields(T,E_t,P_t,g_per,isplot,issmooth)
    %Separate real/imag parts
    if issmooth
        real_E = real(E_t);
        imag_E = imag(E_t);
        real_P = real(P_t);
        imag_P = imag(P_t);

        %smooth real E
        pv_cutoff = 100;    %peak-valley cutoff
        [peaks,valleys] = find_contiguous_regions(real_E,pv_cutoff);
        real_E_sg = filter_selectively(T,real_E,peaks,valleys,isplot);

        %smooth imag E
        [peaks,valleys] = find_contiguous_regions(imag_E,pv_cutoff);
        imag_E_sg = filter_selectively(T,imag_E,peaks,valleys,isplot);

        %smooth real P
        [peaks,valleys] = find_contiguous_regions(real_P,pv_cutoff);
        real_P_sg = filter_selectively(T,real_P,peaks,valleys,isplot);

        %smooth imag P
        [peaks,valleys] = find_contiguous_regions(imag_P,pv_cutoff);
        imag_P_sg = filter_selectively(T,imag_P,peaks,valleys,isplot);

        %Combine real/imag parts
        E_t = complex(real_E_sg,imag_E_sg);
        P_t = complex(real_P_sg,imag_P_sg);
    end
    
    %Get P_dot and then calculate D_t + SG smoothing
    [T,P_dot_sg,P_t,E_t] = get_P_dot_t(T,P_t,E_t);
    D_t = get_D_t(P_dot_sg,P_t,E_t,g_per);
    
    if issmooth
        D_t_sg = sgolayfilt(D_t,2,21);
    else
        D_t_sg = D_t;
    end
    
%     D_t_sg = real(D_t_sg);
%     D_t = real(D_t);
end