#!/bin/bash 
#SBATCH -J hetdex                # Job name
#SBATCH -N 1                     # Total number of nodes requested
##SBATCH -n 1                    # Total number of Tasks
#SBATCH --ntasks 48              # using --ntasks for consistency with spr, but same as --ntasks_per_node
#SBATCH -p skx               # Queue name
#SBATCH -o HETDEX.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -t 13:00:00               # Run time (hh:mm:ss)
#SBATCH -A Hobby-Eberly-Telesco
#SBATCH --mail-user dustin@astro.as.utexas.edu
#SBATCH --mail-type all


module load pylauncher

##export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_WORKDIR=$(pwd)
WORKDIR=$LAUNCHER_WORKDIR
##export LAUNCHER_SCHED=interleaved

module load cdtools
distribute.bash raw/changedate/virus/virus0000changeshot.tar
distribute.bash datafiles/datevshot.dithall
distribute.bash datafiles/datevshotsedtp_f.dat
distribute.bash datafiles/fluxcor.use
distribute.bash datafiles/fwhm.use
distribute.bash datafiles/norm.dat
wait

##export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
##export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE=rt1c
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
##export LAUNCHER_SCHED=interleaved
cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"
run_pylauncher

wait
echo "FIRST PART"

export LAUNCHER_JOB_FILE=rt1
CONTROL_FILE=$LAUNCHER_JOB_FILE
run_pylauncher
wait

runpost
