# /Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nosplash -nodesktop -r "run('blade_runner.m');exit;"

RUN_CMD="/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nosplash -nodesktop -r run('~/Dropbox/dev/cftd-solver/modules/core/blade_runner.m');exit;"
export $RUN_CMD

PROGRAM_NAME="/Applications/MATLAB_R2017a.app/bin/matlab"
export PROGRAM_NAME

MATLAB_SETTINGS=
export MATLAB_SETTINGS

PROGRAM_ARGS="-nodisplay -nosplash -nodesktop -r $(MATLAB_SETTINGS);run('~/Dropbox/dev/cftd-solver/modules/core/blade_runner.m');exit;"
export PROGRAM_ARGS

RUN_CMD="$PROGRAM_NAME $PROGRAM_ARGS"
echo RUN_CMD

PERF_LOG_FILE="outside_log_file.txt"
export PERF_LOG_FILE

MAT_LOG_FILE="inside_mem_log.txt"
export MAT_LOG_FILE

PERF_OUTPUT="cftd-perf.csv"
export PERF_OUTPUT

# Cleanup files
rm -f $PERF_LOG_FILE

echo $PERF_LOG_FILE

interval=1

main() {
    ($RUN_CMD) & pid=$!

    while is_running $pid; do
        report_metrics $pid $PERF_LOG_FILE
        sleep $interval
    done
}

is_running() { (kill -0 ${1:?is_running: missing process ID}) 2>& -; }


function report_metrics {
        TDSS_PID=$1
        TDSS_LOG=$2
        
#       echo $TDSS_PID
        top -l 2 | awk -v ptime="$(date '+%s')" -v tdss="$TDSS_PID" -v logtime="$(date '+%Y-%m-%d,%H:%M:%S')" -v logfile="$TDSS_LOG" -v host=`hostname`  '
              $1 ~ /^[0-9]+$/ && tolower($2) ~ /matlab/ {printf "%s,%s,%s,%s,%s,%s,%s\n", host, $1, logtime, ptime, $3, $8, $2 >> logfile}'
}

main

join -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.4,2.5,2.6,2.8,2.7 -t ','  \
	<(<$PERF_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) \
	<(<$MAT_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) | \
awk -F, '$5 >= $9 && $5 <= $10 {print $0}' > $PERF_OUTPUT


