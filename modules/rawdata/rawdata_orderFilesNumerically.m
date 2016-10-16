function [ordered_files,fnum] = rawdata_orderFilesNumerically(files)
    fnum = zeros(1,length(files));
    for i = 1:length(files)
        name_this = files(i).name;
        num_this  = name_this((find(name_this == '_',1,'last')+1):(find(name_this == '.')-1));
        fnum(i)   = str2double(num_this);
    end
    [fnum,m] = sort(fnum);
    ordered_files = files(m);
end