function param_vecs = setup_genParameterVectors(C_parameters,order)
%Each parameter should be supplied in a separate column of each cell
    R = 1;
    C = 0;
    for i = 1:length(C_parameters)
        [r,c] = size(C_parameters{i});
        C = C + c;
        R = R * r;
    end
    
    param_vecs = zeros(R,C);
    cstart = 1;
    k_prev = R;
    for i = 1:length(C_parameters)
        this_param = C_parameters{i};
        r = size(this_param,1);
        c = size(this_param,2);
        
        k  = k_prev/r;
        kk = R/k_prev;
        
        r_offset = (0:kk-1)*k_prev;
        for ii = 1:length(r_offset)
            rstart = r_offset(ii) + (0:r-1)*k + 1;
            for iii = 1:length(rstart)
                assign_this = repmat(this_param(iii,:),k,1);
                last = rstart(iii) + size(assign_this,1) - 1;
                param_vecs(rstart(iii):last,cstart:(cstart+c-1)) = assign_this;
            end
        end
        cstart = cstart + c;
        k_prev = k;
    end
    param_vecs = param_vecs(:,order);
end
