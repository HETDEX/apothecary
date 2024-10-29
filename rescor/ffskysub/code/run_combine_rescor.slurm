#!/bin/bash 
#SBATCH -J combine                # Job name
#SBATCH -N 1                   # Total number of nodes requested
#SBATCH --ntasks-per-node 1     # Tasks per node
#SBATCH -p development            # Queue name
#SBATCH -o /work/05865/maja_n/ls6/ffskysub/slurm_out/HETDEX.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -t 02:00:00              # Run time (hh:mm:ss)
#SBATCH -A AST23008




eval "$(/work/05865/maja_n/stampede2/anaconda3/bin/conda shell.bash hook); unset PYTHONPATH"
conda activate nbodykit-env
which python


cd /tmp
mkdir maja_n/
mkdir maja_n/ffskysub/


for i in {202201..202212}
do
    python /work/05865/maja_n/ls6/ffskysub/code/combine_rescor_by_month.py $i 1 0 &
done 
wait


rsync -av maja_n/ffskysub/ /scratch/05865/maja_n/ffskysub/

echo "Done copying files."

#python combine_rescor_by_month.py 2018 1 0

#python plot_shot.py 20230826016

echo "Job finished"