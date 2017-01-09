function time_vec = setup_divideTimeVector(t0,t1,cratio)
    time_vec = t0:(cratio*(t1-t0)):t1;
    if time_vec(end) ~= t1
        time_vec = [time_vec t1];
    end
end