function diag_save1DCalcData(S_coredata,data_id,data_type,pstep,issaveCoeffs,isplotCoeffs,issaveField)
    tvec   = S_coredata.tvec;
    sratio = S_coredata.sratio;
    
    [~,CFvecs,~,x]    = cavity_loadBasis(S_coredata.basis_loc,S_coredata.k_a);
    
    save_group0 = diag_getSaveGroup0(tvec,sratio);

    T = [];
    Y = [];
    f=figure;
    c=lines(S_coredata.nCF);
    for i = 1:(length(tvec)-1)
        [T_this,Y_this] = rawdata_load(S_coredata.temp_dir,data_id,data_type,pstep,i);
        
        %Save field
        if issaveField
            T_next      = T_this;
            Yfield_next = userdata_calcTemporalField(Y_this,CFvecs,x,S_coredata.x0);
            rawdata_saveFieldContinuously(T_next,Yfield_next,S_coredata.data_dir,data_id,'field',pstep,i);
            clear Y_next;
        end
        
        %Plot coeffs
        if isplotCoeffs
            for ii = 1:S_coredata.nCF
                plot(T_this,abs(Y_this(:,ii)),'color',c(ii,:)); hold on;
            end
        end
        
        if i >= save_group0 && issaveCoeffs
            T = [T; T_this];
            Y = [Y; Y_this];
        end
    end
    
    if isplotCoeffs
        title([plotting_getStandardTitle(S_coredata) ' pstep=' num2str(pstep) ' (' data_id ',' data_type ')']);
        xlabel('t');
        lgnd = cell(1,S_coredata.nCF);
        for i = 1:S_coredata.nCF
            lgnd{i} = [lower(data_id) '$_{' num2str(i) '}$(t)'];
        end
        legend(lgnd,'interpreter','latex','location','eastoutside','fontsize',14);
        G = get(gcf,'Position');
        set(gcf,'Position',[G(1) G(2) 2*G(3) 2*G(4)]);
    %     export_fig([S_coredata.diag_dir '/' data_id data_type '_pstep_' num2str(pstep) '.png'],'-m2');
        saveas(gcf,[S_coredata.diag_dir '/' data_id data_type '_pstep_' num2str(pstep) '.png']);
    end
    close(f);
    
    %Save coeffs
    if issaveCoeffs
        rawdata_save(T,Y,S_coredata.data_dir,data_id,data_type,pstep);
    end
end