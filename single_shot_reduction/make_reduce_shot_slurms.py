"""
simple script to take all the selected shots (e.g. of a type or in a date range)
and build appropriately limited parallel_XX.run and .slurm files

User should edit near the top to set the email address for notifications
    and alter the selection as needed

"""

import numpy as np
import os
from astropy.table import Table
import shutil

NOTIFICATION_EMAIL = None # "XXX@astro.as.utexas.edu"
SU_ACCOUNT = "AST23008"
BASE_FILENAME = "parallel"
MIN_TOTAL_EXPTIME = 300.0 #seconds
MAX_TOTAL_EXPTIME = 999999.0 #seconds
START_DATE = 20240801  #inclusive
END_DATE = 20240831    #inclusive

QUEUE = "normal"
MAX_NODES = 5 #this is a TACC request to limit
MAX_TASKS_PER_NODE = 11 #this is a memory issue
NOMINAL_TIME_REQUEST = "3:00:00"  #this per SLURM file, the assumption here is there are 5*11 = 55 shots per file
#for future, rather than fix a time, let it be estimated, like in elixer
#PER_REDUCTION_TIME = 3.0 #hours (nominal worst case ... a 3-dither, deep observation)
NUM_SLURM_FILES = None # if None, let the script fit, otherwise use this value?

#if you have defined a common shell wrapper, use it, otherwise use the direct python call to the python script
if shutil.which("reduce_shot"):
    BASE_CMD = "reduce_shot "
else:
    BASE_CMD = "python reduce_shot.py "

EXTRA_SWITCHES = "--linedet_filter 0 " #"--multifits_only --clean 0"
#EXTRA_SWITCHES = "--multifits_only --local_het_raw_path /scratch/03261/polonius/het_raw --clean 0"

#note: these are the defaults
#--cal_flux 'sdss' --cal_astro 'gaia'  --clean 1


#read the table
if os.path.exists("virus_shot_summary_table.tab"):
    T = Table.read("virus_shot_summary_table.tab",format="ascii")
elif os.path.exists("/work/03261/polonius/hetdex/single_shot/virus_shot_summary_table.tab"):
    T = Table.read("/work/03261/polonius/hetdex/single_shot/virus_shot_summary_table.tab", format="ascii")
else:
    print("Fatal. Cannot locate virus_shot_summary_table.tab")
    exit(-1)


##################
# shot selection
##################
sel = np.array(T['exptime'] >= MIN_TOTAL_EXPTIME) * np.array(T['exptime'] <= MAX_TOTAL_EXPTIME)

#new version of the virus_shot_summary_table.tab is a whitelist, so no need to kick out VDFI as it is not already there
#still, just to be safe ....
#sel = sel * np.array(T['name'] == 'parallel')
sel = sel * np.array([n[0:4] != 'VDFI' for n in T['name']])
sel = sel * np.array(T['date'] >= START_DATE) * np.array(T['date'] < END_DATE)
#sel = sel * np.array(T['num_exp'] == 1)
#T[sel]

def find_min_product(max_a: int, max_b: int, c: int):
    """
    Find integers a <= max_a and b <= max_b such that a * b is the smallest
    value >= c.

    Returns (a, b), or raises ValueError if no valid pair exists.
    """
    if max_a * max_b < c:
        return max_a, max_b

    best_a, best_b, best_product = None, None, float("inf")

    for a in range(1, max_a + 1):
        # Smallest b such that a * b >= c
        b = -(-c // a)          # ceiling division: equivalent to math.ceil(c / a)
        if b < 1:
            b = 1
        if b > max_b:
            continue            # this value of a can't satisfy the constraint
        product = a * b
        if product < best_product:
            best_product, best_a, best_b = product, a, b

    if best_a is None:
        return max_a, max_b

    return best_a, best_b

def make_slurm(basename, max_nodes, max_tpn, tasks):

    global NOTIFICATION_EMAIL, QUEUE, NOMINAL_TIME_REQUEST, SU_ACCOUNT

    if tasks >= (max_nodes * max_tpn):
        tpn = max_tpn
        nodes = max_nodes
    else:
        nodes, tpn = find_min_product(max_nodes,max_tpn,tasks)

    #nodes = tasks // max_tpn + (1 if tasks % max_tpn !=0 else 0)
    print(f"nodes, ntasks_per_node, total tasks: {nodes}, {tpn}, {tasks}")

    #nodes = min(nodes, max_nodes)
    jobfile = f"{basename}.run"
    slurmfile = f"{basename}.slurm"

    slurmstr = f"""#!/bin/bash
#SBATCH -J HETSSR                # Job name
#SBATCH -N {nodes}                     # Total number of nodes requested
#SBATCH --ntasks-per-node {tpn}     # tasks per node
#SBATCH -p {QUEUE}                # Queue name
#SBATCH -o HETSSR.o%j            # Name of stdout output file (%j expands to jobid)
#SBATCH -t {NOMINAL_TIME_REQUEST}               # Run time (hh:mm:ss)
#SBATCH -A {SU_ACCOUNT}
"""
    if NOTIFICATION_EMAIL is not None:
        slurmstr += f"#SBATCH --mail-user {NOTIFICATION_EMAIL} \n"
        slurmstr += "#SBATCH --mail-type all \n"
    slurmstr += f"""module load launcher
export LAUNCHER_PLUGIN_DIR=$TACC_LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_WORKDIR=$(pwd)
export LAUNCHER_JOB_FILE={jobfile}
WORKDIR=$LAUNCHER_WORKDIR
CONTROL_FILE=$LAUNCHER_JOB_FILE
export LAUNCHER_SCHED=interleaved

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun

echo " "
echo " Parameteric Job Complete"
echo " "
    """

    with open(slurmfile,"w") as f:
        f.write(slurmstr)
        f.write("\n")
    print(f"Wrote: {slurmfile}")


############
# main
############

nodes = MAX_NODES
tasks_per_node = MAX_TASKS_PER_NODE
if "multifits_only" in EXTRA_SWITCHES:
    step = len(T[sel])
    num_files = 1
elif NUM_SLURM_FILES is not None:
    if NUM_SLURM_FILES > len(T[sel]):
        step = 1
        num_files = len(T[sel])
    else:
        step = len(T[sel]) // NUM_SLURM_FILES
        num_files = len(T[sel]) // step + (1 if len(T[sel]) % step != 0 else 0)
else:
    step = nodes * tasks_per_node  # how many per file, 55 is 5 nodes with 11 per node (given the memory requirements)
    num_files = len(T[sel]) // step + (1 if len(T[sel]) % step != 0 else 0)

print(f"Making {num_files} files at {tasks_per_node} tasks_per_node and {nodes} nodes")

start_ix = -step
stop_ix = 0
for fi in range(num_files):
    start_ix += step
    stop_ix = min(start_ix + step, len(T[sel]))
    rn = str(fi + 1).zfill(2)
    num_lines = stop_ix - start_ix
    if NOTIFICATION_EMAIL is not None:
        lines = [f"{BASE_CMD} --email {NOTIFICATION_EMAIL} {s}\n" for s in T['shotid'][sel][start_ix:stop_ix]]
    else:
        lines = [f"{BASE_CMD} {EXTRA_SWITCHES} {s}\n" for s in T['shotid'][sel][start_ix:stop_ix]]

    basename = f"{BASE_FILENAME}_{rn}"
    fn = f"{basename}.run"
    with open(fn, "w") as f:
        for i, row in enumerate(lines):
            f.write(lines[i])

        print(f"Wrote: {fn} ; {len(lines)} lines")

    make_slurm(basename, nodes, tasks_per_node, num_lines)
