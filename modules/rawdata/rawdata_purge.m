function rawdata_purge
    datafiles = dir('num*');
    for i = 1:length(datafiles)
        rmdir(datafiles(i).name,'s')
    end
    if exist('results','dir') == 7
        rmdir('results','s');
    end
    if exist('integrals','dir') == 7
        rmdir('integrals','s');
    end
    rmdir('slurms','s')
    delete for_bash
    delete TD_init_parameters.mat
    delete TD_D0_parameters.m~
end