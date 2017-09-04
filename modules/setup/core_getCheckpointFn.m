function checkpoint_fn =  core_getCheckpointFn(pstep,cp_dir)
    checkpoint_fn = [cp_dir '/checkpoint_p' num2str(pstep) '.mat'];
end