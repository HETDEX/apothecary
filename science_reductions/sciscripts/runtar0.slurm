#!/bin/bash
#
#------------------Scheduler Options--------------------
#SBATCH -J HETDEX              # Job name
#SBATCH -n 1                   # Total number of tasks
#SBATCH -N 1                   # Total number of Nodes
#SBATCH -p hadoop              # Queue name
#SBATCH -o HETDEX.o%j          # Name of stdout output file (%j expands to jobid)
#SBATCH -t 09:00:00            # Run time (hh:mm:ss)
#SBATCH -A Hobby-Eberly-Telesco
#------------------General Options---------------------
module use /work/01255/siliu/launcher/modulefiles/
module load launcher/3.6
module use /work/01255/siliu/stampede2/ooops/modulefiles/
module use /work/00410/huang/share/patch/modules
module unload xalt

echo $SLURM_JOB_NODELIST
hostname
pwd

export TACC_LAUNCHER_PPN=24
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=$(pwd)
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_SCHED=interleaved
export LAUNCHER_JOB_FILE=runtar0

#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun SLURM $EXECUTABLE $WORKDIR $CONTROL_FILE $PHI_WORKDIR $PHI_CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
