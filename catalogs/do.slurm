#!/bin/bash
#SBATCH -J HETCAT                # Job name
#SBATCH -N 1                     # Total number of nodes requested
#SBATCH --ntasks-per-node 50     # tasks per node
#SBATCH -p normal                # Queue name
#SBATCH -o HETCAT.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -e HETCAT.e%j            # Name of stderr output file (%j expands to jobid)
#SBATCH -t 5:00:00               # Run time (hh:mm:ss)
#SBATCH -A AST23008
#SBATCH --mail-user dustin@astro.as.utexas.edu
#SBATCH --mail-type all


module load launcher
export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE=do.run
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
export LAUNCHER_SCHED=interleaved

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun

echo " "
echo " Parameteric Job Complete"
echo " "
