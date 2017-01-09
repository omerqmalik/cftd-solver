function save_group0 = diag_getSaveGroup0(tvec,sratio)
    save_t0 = round(tvec(end)*(1 - sratio));
    tvec_new = sort([tvec save_t0]);
    pos_new = find(tvec_new == save_t0);
    if length(pos_new) > 1
        save_group0 = pos_new(1);
    else
        if pos_new <= 2
            save_group0 = 1;
        elseif pos_new >= length(tvec)
            save_group0 = length(tvec) - 1;
        else
            save_group0 = pos_new - 1;
        end
    end
end