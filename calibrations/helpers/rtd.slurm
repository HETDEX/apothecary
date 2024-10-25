#!/bin/bash 
#SBATCH -J HETDEX                # Job name
#SBATCH -N 1                     # Total number of nodes requested
#SBATCH -n 1                     # Total number of Tasks
#SBATCH -p development           # Queue name
#SBATCH -o HETDEX.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -t 02:00:00              # Run time (hh:mm:ss)
#SBATCH -A AST23008
module load launcher
export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE=rt
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
export LAUNCHER_SCHED=interleaved

#----------------
# Job Submission
#----------------

cd $WORKDIR/ 
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun

echo " "
echo " Parameteric Job Complete"
echo " " 
