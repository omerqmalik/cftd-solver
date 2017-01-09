function S_setupdata = setup_loadParameters(cav_dir,num)
    load([cav_dir '/TD_init_parameters'],'S_setupdata');
    load([cav_dir '/TD_init_parameters'],'S_pumpdata');
    if nargin > 1
        S_setupdata = S_setupdata(num);
    end
    [S_setupdata.pump_type] = deal(S_pumpdata.type);
end