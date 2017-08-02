setenv ( 'CFTD_TEMP_PATH', 'C:\Users\nvadmin\Dropbox\workspace\cftd-solver');

TD_D0_parameters;

mex -I/Users/tanay/Dropbox/dev/boost_1_64_0/ runge_kutta4.cpp;

clearvars -global
set_defaults(); %sets the values of global variables
global ode_options;
ode_options.use_ode_cpp = true;

core_runTDSS('.',1,1);

