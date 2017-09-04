function diag_save2DCalcData(S_coredata,data_id,data_type,pstep)
    tvec   = S_coredata.tvec;
    sratio = S_coredata.sratio;
    
    save_group0 = diag_getSaveGroup0(tvec,sratio);

    Y = zeros(S_coredata.nCF^2,1);
    T = [];
    for i = save_group0:(length(tvec)-1)
        [T_this,Y_this] = rawdata_load(S_coredata.temp_dir,data_id,data_type,pstep,i);
        T = [T T_this];
        Y = Y + Y_this;
    end
    T = sort(unique(T));
    T = [T(1) T(end)];
    Y = Y/length(save_group0:(length(tvec)-1));
    Y = reshape(Y,[S_coredata.nCF,S_coredata.nCF]).';
    
    f = plotting_makeFig;
    imagesc(Y);
    colormap hot;
    title([plotting_getStandardTitle(S_coredata) ' pstep=' num2str(pstep) ' (' data_id ',' data_type ')']);
    
    h1 = colorbar;
    ylabel(h1,['$\bar{D}_{mn}(t)$'],'interpreter','latex');
    
    saveas(f,[S_coredata.diag_dir '/' data_id data_type '_' num2str(pstep) '.png']);
    
    userdata_saveT(T,S_coredata.data_dir,data_id,data_type);
    userdata_saveY(Y,S_coredata.data_dir,data_id,data_type,pstep);
    close(f);
end