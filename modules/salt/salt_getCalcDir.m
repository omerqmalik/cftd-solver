function calc_dir = salt_getCalcDir(cav_dir,num,prel_f,sz,S_setupdata)
    calc_dir = [cav_dir '/' S_setupdata(num).calc_dir '/salt/prel_1-' num2str(prel_f) '_sz_' num2str(sz)];
end