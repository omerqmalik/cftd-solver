function set_defaults()
    global ode_options;
    global mem_log_file;
    
    mem_log_file = "/tmp/inside_mem_log.txt";
    
    %Clear out mem_log_file
    temp_fid = fopen(mem_log_file, 'wt');
    fclose(temp_fid);
    
    
    % log level are "ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"
    ode_options.log_level = 'INFO';
    ode_options.use_ode_cpp = false;
    ode_options.odeint_const = true;
    ode_options.dt = 0.1;
    ode_options.abs_error = 1.0e-6;
    ode_options.rel_error = 1.0e-3;
    
    global benchmarking;
    benchmarking = true;
end


