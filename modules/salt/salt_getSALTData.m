function [k,a,a2,pump] = salt_getSALTData(cav_dir,num,prel_f,sz)
    S_coredata = salt_loadParameters(cav_dir,num,prel_f,sz);
    
    datafiles = dir([S_coredata.data_dir '/nlsolver_*.mat']);
    datafiles = rawdata_orderFilesNumerically(datafiles);
    len_DF = length(datafiles);
    
    k = zeros(len_DF,S_coredata.nCF);
    a = zeros(S_coredata.nCF,S_coredata.nCF,len_DF);
    pump = zeros(1,len_DF);
    
    for i = 1:len_DF
        load([S_coredata.data_dir '/' datafiles(i).name]);
        
        if i == 1
            k(1,1)   = S_nlsol.k;
            a(:,1,1) = S_nlsol.a;
            
            kprev = S_nlsol.k;
        else
            ind = 1:length(S_nlsol.k);
            for ii = 1:length(kprev)
                [~,minpos] = min(abs(kprev(ii)-S_nlsol.k));
                k(i,ii) = S_nlsol.k(minpos);
                a(:,ii,i) = S_nlsol.a(:,minpos);
                
                ind(ind == minpos) = [];
            end
            if length(S_nlsol.k) > length(kprev)
                k(i,length(S_nlsol.k)) = S_nlsol.k(ind);
                a(:,length(S_nlsol.k),i) = S_nlsol.a(:,ind);
            end
            kprev = k(i,k(i,:)>0);
        end
        pump(i) = S_nlsol.d0_this;
    end
    
    ind = sum(k,1) > 0;
    k = k(:,ind);
    a = a(:,ind,:);
    
    a2 = transpose(squeeze(sum(abs(a).^2,1)));
end