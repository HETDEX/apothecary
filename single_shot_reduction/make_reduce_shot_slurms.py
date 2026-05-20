"""
simple script to take all the selected shots (e.g. of a type or in a date range)
and build appropriately limited parallel_XX.run and .slurm files

User should edit near the top to set the email address for notifications
    and alter the selection as needed

"""

import numpy as np
from astropy.table import Table

NOTIFICATION_EMAIL = "dustin@astro.as.utexas.edu"
MAX_NODES = 5 #this is a TACC request to limit
MAX_TASKS_PER_NODE = 11 #this is a memory issue
QUEUE = "normal"
NOMINAL_TIME_REQUEST = "4:00:00"
SU_ACCOUNT = "AST25019"
BASE_CMD = "python reduce_shot.py --cal_flux 'sdss' --cal_astro 'gaia' --clean 1"
BASE_FILENAME = "parallel"

#read the table
T = Table.read("virus_shot_summary_table.tab",format="ascii")

##################
# shot selection
##################
sel = np.array(T['exptime'] > 300) * np.array(T['exptime'] < 9999)

sel = sel * np.array(T['name'] == 'parallel')
sel = sel * np.array([n[0:4] != 'VDFI' for n in T['name']])
sel = sel * np.array(T['date'] >= 20260101) * np.array(T['date'] < 20260201)
#sel = sel * np.array(T['num_exp'] == 1)
#T[sel]


def make_slurm(basename, max_nodes, max_tpn, tasks):

    global NOTIFICATION_EMAIL, QUEUE, NOMINAL_TIME_REQUEST, SU_ACCOUNT

    tpn = tasks if tasks < max_tpn else max_tpn

    nodes = tasks // max_tpn + (1 if tasks % max_tpn !=0 else 0)
    print(nodes, tasks, max_tpn)

    nodes = min(nodes, max_nodes)
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
#SBATCH --mail-user {NOTIFICATION_EMAIL}
#SBATCH --mail-type all


module load launcher
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
step = nodes * tasks_per_node  # how many per file, 55 is 5 nodes with 11 per node (given the memory requirements)
num_files = len(T[sel]) // step + 1 if len(T[sel]) % step != 0 else 0

print(f"Making {num_files} files at {tasks_per_node} tasks_per_node and {nodes} nodes")

start_ix = -step
stop_ix = 0
for fi in range(num_files):
    start_ix += step
    stop_ix = min(start_ix + step, len(T[sel]))
    rn = str(fi + 1).zfill(2)
    num_lines = stop_ix - start_ix
    lines = [f"{BASE_CMD} --email {NOTIFICATION_EMAIL} {s}\n" for s in T['shotid'][sel][start_ix:stop_ix]]

    basename = f"{BASE_FILENAME}_{rn}"
    fn = f"{basename}.run"
    with open(fn, "w") as f:
        for i, row in enumerate(lines):
            f.write(lines[i])

        print(f"Wrote: {fn} ; {len(lines)} lines")

    make_slurm(basename, nodes, tasks_per_node, num_lines)