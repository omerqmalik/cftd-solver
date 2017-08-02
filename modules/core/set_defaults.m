function set_defaults()
    global ode_options;
    ode_options.log_level = 'Debug';
    ode_options.use_ode_cpp = true;
    ode_options.odeint_const = true;
    
    global benchmarking;
    benchmarking = true;
end


