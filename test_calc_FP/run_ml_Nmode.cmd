#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem-per-cpu=56000
#SBATCH --array=0-29
#SBATCH -o "slurms/slurm-%A_%a.out"

cavloc=$(pwd)
pgroups=$(cat "$cavloc/for_bash")
cnum=$[SLURM_ARRAY_TASK_ID/pgroups+1]
pnum=$[SLURM_ARRAY_TASK_ID%pgroups+1]

/usr/licensed/bin/matlab -nodisplay -nosplash -r "addpath '/tigress/omalik/Time Dynamics/cftd-solver/modules/core'; setenv('CFTD_TEMP_PATH','/scratch/gpfs/omalik'); core_runTDSS('$cavloc',$cnum,$pnum)"
