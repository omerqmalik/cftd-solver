function [pump_ind,pumpgrp_ind] = pump_getPumpForTDSSolver(calc_times,pump_D0,pumpgrp_len,pump_type,pumpgrp_num)
    last_pstep = sum(sum(calc_times,1) > 0);
    if last_pstep == size(calc_times,2)
        pump_ind = [];
        pumpgrp_ind = [];
        return;
    end
    
    pumpgrp_start = pumpgrp_len*(pumpgrp_num-1) + 1;
    pumpgrp_end   = pumpgrp_start + pumpgrp_len - 1;
    if pumpgrp_end > length(pump_D0)
        pumpgrp_end = length(pump_D0);
    end
    
    if strcmp(pump_type,'simple')
        pump_ind    = (pumpgrp_start + last_pstep):pumpgrp_end;
        pumpgrp_ind = (last_pstep+1):size(calc_times,2);
    elseif strcmp(pump_type,'hysteresis')
        pump_ind    = (pumpgrp_start + last_pstep):length(pump_D0);
        pumpgrp_ind = pump_ind;
    end
end