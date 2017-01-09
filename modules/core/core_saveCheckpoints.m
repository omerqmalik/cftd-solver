function core_saveCheckpoints(t,checkpoints,checkpoints_fn)
    [t_prev, cp_prev] = core_loadCheckpoints(checkpoints_fn);
    checkpoints = [cp_prev checkpoints];
    t = [t_prev t];
    save(checkpoints_fn,'t','checkpoints');
end