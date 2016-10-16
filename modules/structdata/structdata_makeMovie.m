function S_movie = structdata_makeMovie(movie_type,S_dataArray,S_figArray)
    if strcmp(movie_type,'fft')
        S_movie.type = 'fftmov';
        S_movie.islin = S_figArray(1).islin;
        S_movie.isdec = S_figArray(1).isdec;
        S_movie.fname = S_dataArray(1).fname;
        S_movie.M = plotting_makeMovie(@(S_thisdata,S_thisfig) structdata_plotFFT(S_thisdata,S_thisfig),S_dataArray,S_figArray);
    elseif strcmp(movie_type,'tw')
        S_movie.type = 'twmov';
        S_movie.fname = S_dataArray(1).fname;
        S_movie.M = plotting_makeMovie(@(S_thisdata,S_thisfig) structdata_plotTemporalWaveform(S_thisdata,S_thisfig),S_dataArray,S_figArray);
    end
end