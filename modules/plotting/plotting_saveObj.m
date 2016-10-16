function plotting_saveObj(S_obj,results_dir)
    if strcmp(S_obj.type,'3Dfft') || strcmp(S_obj.type,'tw')
        savefig(S_obj.f,[results_dir '/' S_obj.type '_' S_obj.fname '.fig']);
    elseif strcmp(S_obj.type,'2Dfft')
        results_dir = [results_dir '/' S_obj.type '_' S_obj.fname];
        if ~S_obj.islin
            results_dir = [results_dir '_log'];
        end
        if S_obj.isdec
            results_dir = [results_dir '_d'];
        end
        savefig(S_obj.f,[results_dir '.fig']);
    elseif strcmp(S_obj.type,'fftmov')
        M = S_obj.M;
        results_dir = [results_dir '/' S_obj.type '_' S_obj.fname];
        if ~S_obj.islin
            results_dir = [results_dir '_log'];
        end
        if S_obj.isdec
            results_dir = [results_dir '_d'];
        end
        save([results_dir '.mat'],'M');
    elseif strcmp(S_obj.type,'twmov')
        M = S_obj.M;
        results_dir = [results_dir '/' S_obj.type '_' S_obj.fname];
        save([results_dir '.mat'],'M');
    end
end