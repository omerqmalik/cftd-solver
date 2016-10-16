function S_setupdata = setup_loadParameters(cav_dir,num)
    load([cav_dir '/TD_init_parameters'],'S_setupdata');
    if nargin > 1
        S_setupdata = S_setupdata(num);
    end
end