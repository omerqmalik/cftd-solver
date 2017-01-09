function [t,checkpoints] = core_loadCheckpoints(checkpoints_fn)
    if exist(checkpoints_fn,'file') == 2
        load(checkpoints_fn,'t','checkpoints');
    else
        t = [];
        checkpoints = [];
    end
end