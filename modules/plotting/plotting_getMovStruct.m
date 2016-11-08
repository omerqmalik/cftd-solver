function S_movie = plotting_getMovStruct(movie_type,S_dataArray,opt_islin,opt_isdec)
    S_movie.type = movie_type;
    if strcmp(movie_type,'2DfieldFFTmov')
        S_movie.islin = opt_islin;
        S_movie.isdec = opt_isdec;
    end
    
    fname = [movie_type '_' S_dataArray(1).id S_dataArray(1).type '_p' num2str(S_dataArray(1).pump/S_dataArray(1).th) '-' num2str(S_dataArray(end).pump/S_dataArray(end).th) '_t' num2str(S_dataArray(1).t(1)) '-' num2str(S_dataArray(1).t(end))];
    
    append_str = '';
    if nargin > 2
        if ~opt_islin
            append_str = [append_str '_log'];
        end
        if opt_isdec
            append_str = [append_str '_d'];
        end
    end
    S_movie.fname = [fname append_str];
end