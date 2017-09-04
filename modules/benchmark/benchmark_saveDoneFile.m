function benchmark_saveDoneFile(saveloc,pstep)
    filename = [saveloc '/done_' num2str(pstep)];
    done = 1;
    save(filename,'done');
end