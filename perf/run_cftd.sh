#
# RUN Script for CFTD Solver
# All parameters should be updated in run_cftd_params.sh file
#

#source 
. ./run_cftd_params.sh

RUN_CMD="$PROGRAM_NAME $PROGRAM_ARGS"
echo RUN_CMD

# Cleanup files
PERF_LOG_FILE="$LOG_DIR $METRICS_LOG_FILE"
echo $PERF_LOG_FILE
rm -f $PERF_LOG_FILE

MAT_LOG_FILE="$LOG_DIR $MATLAB_LOG_FILE"
echo $MAT_LOG_FILE

PERF_OUTPUT="$LOG_DIR $PERF_OUTPUT"

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
        TDSS_PID=$1 # unused for now
        TDSS_LOG=$2
        
        top -l 2 | awk -v ptime="$(date '+%s')" -v logtime="$(date '+%Y-%m-%d,%H:%M:%S')" -v logfile="$TDSS_LOG" -v host=`hostname`  '
              $1 ~ /^[0-9]+$/ && tolower($2) ~ /matlab/ {printf "%s,%s,%s,%s,%s,%s,%s\n", host, $1, logtime, ptime, $3, $8, $2 >> logfile}'
}

main

join -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.4,2.5,2.6,2.8,2.7 -t ','  \
	<(<$PERF_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) \
	<(<$MAT_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) | \
awk -F, '$5 >= $9 && $5 <= $10 {print $0}' > $PERF_OUTPUT


