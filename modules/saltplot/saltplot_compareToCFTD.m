function saltplot_compareToCFTD(cav_dir,num,prel_f,sz,beta)
    [k,a,a2,pump] = salt_getSALTData(cav_dir,num,prel_f,sz);
    pump = pump/min(pump);
    
    if nargin > 4
        [~,ind] = helpers_getClosestMatch(beta,pump);
        k = k(ind,:);
        if size(a2,1) == 1
            a2 = a2.';
        end
        a2 = a2(ind,:);

        userplot_plotTWCoeffs(cav_dir,num,beta);
        xL = get(gca,'xlim');

        c = colormap(lines(length(k)));

        for i = 1:length(k)
            line(xL,[a2(i) a2(i)],'color',c(i,:),'linestyle','--');
        end
    else
        S_EcoeffsAll   = structdata_load(cav_dir,num,'E','coeffs',1,0);
        
        Y_a = squeeze(mean(abs(S_EcoeffsAll.Y).^2,1));
        if size(Y_a,2) == 1
            Y_a = Y_a.';
        end
        Y_t = S_EcoeffsAll.Y;
        
        CFTD_pump  = S_EcoeffsAll.pump/S_EcoeffsAll.th;
        CFTD_pumpt = CFTD_pump(CFTD_pump <= prel_f);
        Y_a = Y_a(:,CFTD_pump <= prel_f);
        Y_t = Y_t(:,:,CFTD_pump <= prel_f);
        
        figure;
        plot(CFTD_pumpt,Y_a); hold on;
        set(gca,'ColorOrderIndex',1);
        plot(pump,a2,'linestyle','none','marker','o');
        
%         %%%% Compare to analytic %%%%
%         load TD_init_parameters.mat
%         load integrals/TD_integrals_nCF_1_ka_20.9433_n_3.mat
%         A = real(A);
%         om = -S_setupdata.g_per*imag(S_setupdata.CFvals)^2/(real(S_setupdata.CFvals)*imag(S_setupdata.CFvals) + S_setupdata.k_a*S_setupdata.g_per)/2;
%         gc = S_setupdata.g_per^2/(S_setupdata.g_per^2 + om^2);
%         I = (S_setupdata.pump/S_setupdata.th - 1)/A/gc;
%         plot(S_setupdata.pump/S_setupdata.th,I,'o');
%         %%%%%%%
        
        figure;
        [w,~,fft_w_mag] = userdata_calcFFT(S_EcoeffsAll.t,Y_t,S_EcoeffsAll.rframe);
        [~,max_ind] = max(fft_w_mag.^2,[],1);
        max_locs = w(squeeze(max_ind));
        plot(CFTD_pumpt(2:end),max_locs(:,2:end)); hold on;
        set(gca,'ColorOrderIndex',1);
        
        for i = 1:size(k,2)
            plot(pump(k(:,i)>0),k(k(:,i)>0,i),'linestyle','none','marker','o');
        end
    end
end