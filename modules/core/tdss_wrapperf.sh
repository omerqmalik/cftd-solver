program_name="$1"
args_for_prog=
log_file=$2

max=30 interval=1 n=0

main() {
    ($program_name $args_for_prog) & pid=$!

    while is_running $pid; do
        report_metrics $pid $log_file
        sleep $interval
    done
#    echo
}

is_running() { (kill -0 ${1:?is_running: missing process ID}) 2>& -; }


function report_metrics {
        TDSS_PID=$1
        TDSS_LOG=$2
        echo $TDSS_PID
        top -l 2 | awk -v tdss="$TDSS_PID" -v logtime="$(date +%T)" -v logfile="$TDSS_LOG" '
                $1 == tdss {printf "%s,%s,%s,%s,%s\n", logtime, $1, $3, $8, $2 >> logfile}'
}

main