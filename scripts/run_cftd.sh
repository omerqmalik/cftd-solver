#
# RUN Script for CFTD Solver on OS X
# All parameters should be updated in run_cftd_params.sh file
#

# source 
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
echo $PERF_OUTPUT

# Sleep duration
interval=1

main() {
    ($RUN_CMD) & pid=$!

    while is_running $pid; do
        record_metrics $pid $PERF_LOG_FILE
        sleep $interval
    done
}

is_running() { (kill -0 ${1:?is_running: missing process ID}) 2>& -; }


function record_metrics {
        TDSS_PID=$1 # unused for now
        TDSS_LOG=$2
        

        top -l 2 | awk -v ptime="$(date '+%s')" -v logtime="$(date '+%Y-%m-%d,%H:%M:%S')" -v logfile="$TDSS_LOG" -v host=`hostname`  '
              $1 ~ /^[0-9]+$/ && tolower($2) ~ /matlab/ {printf "%s,%s,%s,%s,%s,%s,%s\n", host, $1, logtime, ptime, $3, $8, $2 >> logfile}'

        # top -l 2 displays cpu and memory stats for every process
        # awk -v ptime: number of seconds from epoch
        #     -v logtime: entire date to the second
        #     -v logfile: log file from input to funtion
        #     -v host: hostname retrieved from system command
        #
        #     $1 ~ /^[0-9]+$/: checks if first term is a number
        #                      Specifically, a process ID
        #
        #     tolower($2) ~ /matlab/: checks if the process name
        #                             has "matlab" in it
        #     printf: print fields and variables to logfile
        #               variables: 
        #                   host, logtime, ptime
        #               fields:
        #                   $1: process ID
        #                   $3: %CPU 
        #                   $8: memory
        #                   $2: command name
        #

}

main

join -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.4,2.5,2.6,2.8,2.7 -t ','  \
	<(<$PERF_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) \
	<(<$MAT_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1) | \
awk -F, '$5 >= $9 && $5 <= $10 {print $0}' > $PERF_OUTPUT

# join -j1: joins files on the first field
#       -o1.2...: output fields from files 1 and 2
#                 i.e. 1.3 refers to the 3rd field from file 1
#
#       -t ",": delimited by comma
#
#
# (<$PERF_LOG_FILE awk -F, '{print $1"-"$2","$0}' | sort -k1,1)
# This is file 1
#
# awk -F, : parses string using comma delimiter
#
# print $1"-"$2","$0: combines fields 1 and 2 into a new field, used to join
#
# sort -k1,1: sorts the file based on the first field
#             this is required by join
#
#
# File 2 follows same structure as file 1
#
#
# awk -F, '$5 >= $9 && $5 <= $10 {print $0}' > $PERF_OUTPUT
# conditionally writes the output of join to output
#
# awk -F, : using comma delimiter
#
# '$5 >= $9 && $5 <= $10 : the condition to be printed: 
#                          if $5 (the time) is both
#                          greater than the start time ($9)
#                          and less than the end time ($10)
#
# print $0: prints the entire line
#














