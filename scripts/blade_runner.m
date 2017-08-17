setenv ('CFTD_TEMP_PATH', '../');
addpath(genpath('../'));

mex -I../../boost_1_64_0/ ../modules/core/runge_kutta4.cpp; 

TD_D0_parameters;

clearvars -global
set_defaults(); %sets the values of global variables
global ode_options;
ode_options.use_ode_cpp = true;

core_runTDSS('.',1,1);
