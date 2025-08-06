#!/bin/bash 
#SBATCH -J HETDEX                # Job name
##SBATCH -N 3                     # Total number of nodes requested
#SBATCH -N 4                     # Total number of nodes requested
#SBATCH -n 80                    # Total number of Tasks
#SBATCH -p normal               # Queue name
#SBATCH -o HETDEX.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -t 04:00:00               # Run time (hh:mm:ss)
##SBATCH -A AST23008
#SBATCH -A Hobby-Eberly-Telesco
#SBATCH --mail-user dustin@astro.as.utexas.edu
#SBATCH --mail-type all
module load launcher
export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE=template
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
export LAUNCHER_SCHED=interleaved

cd $WORKDIR/ 
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun

if [ -e $WORKDIR/runpost ]
then
   $WORKDIR/runpost
else
   echo "no post processing"
fi

echo " "
echo " Parameteric Job Complete"
echo " " 
