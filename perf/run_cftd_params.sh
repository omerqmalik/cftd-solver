
RUN_CMD="/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nosplash -nodesktop -r run('~/Dropbox/dev/cftd-solver/modules/core/blade_runner.m');exit;"
export RUN_CMD

PROGRAM_NAME="/Applications/MATLAB_R2017a.app/bin/matlab"
export PROGRAM_NAME

#LOG_DIR="~/Dropbox/dev/cftd-solver/perf/"
LOG_DIR=
export LOG_DIR

METRICS_LOG_FILE="outside_log_file2.txt"
export METRICS_LOG_FILE

MATLAB_LOG_FILE="inside_mem_log.txt"
export MATLAB_LOG_FILE

PERF_OUTPUT="cftd-perf2.csv"
export PERF_OUTPUT

#MATLAB_SETTINGS="global mem_log_file; mem_log_file=$MAT_LOG_FILE;"
MATLAB_SETTINGS=
export MATLAB_SETTINGS

PROGRAM_ARGS="-nodisplay -nosplash -nodesktop -r $MATLAB_SETTINGS run('~/Dropbox/dev/cftd-solver/modules/core/blade_runner.m');exit;"
export PROGRAM_ARGS

