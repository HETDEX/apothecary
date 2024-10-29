#!/bin/bash 
#SBATCH -J rescor                # Job name
#SBATCH -N 1                     # Total number of nodes requested
#SBATCH --ntasks-per-node 12     # Tasks per node
#SBATCH -p normal                # Queue name
#SBATCH -o HETDEX_rescor.o%j     # Name of stdout output file (%j expands to jobid)
#SBATCH -t 4:00:00               # Run time (hh:mm:ss)
#SBATCH -A AST23008
#SBATCH --mail-user dustin@astro.as.utexas.edu
#SBATCH --mail-type all

module load launcher
export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE=rescor_combine.run
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
export LAUNCHER_SCHED=interleaved


#----------------
# Error Checking
#----------------

if [ ! -d $WORKDIR ]; then 
        echo " " 
	echo "Error: unable to change to working directory." 
	echo "       $WORKDIR" 
	echo " " 
	echo "Job not submitted."
	exit
fi

if [ ! -x $EXECUTABLE ]; then
   echo " "
   echo "Error: unable to find launcher executable $EXECUTABLE."
   echo " "
   echo "Job not submitted."
   exit
fi

if [ ! -e $CONTROL_FILE ]; then # used to be $WORKDIR/$CONTROL_FILE
   echo " "
   echo "Error: unable to find input control file $CONTROL_FILE."
   echo " "
   echo "Job not submitted."
   exit
fi

#----------------
# Job Submission
#----------------

#mkdir -p /tmp/polonius/ffskysub/improved_spectra/pytables/nomask
#cd /tmp/polonius

$TACC_LAUNCHER_DIR/paramrun

#rsync -a polonius/ffskysub/ /scratch/03261/polonius/red1/rescor

#echo "Done copying files."

echo " "
echo " Parameteric Job Complete"
echo " "
