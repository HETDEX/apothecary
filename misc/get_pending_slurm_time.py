###!/usr/bin/env python
"""
get the total pending SLURM time for an account

(e.g. use what time is still outstanding here + used time from TACC portal to figure what time is left)

"""


import numpy as np
import subprocess


project = "AST23008"

r = subprocess.run(f"squeue -o \"%a %j %u %2t %4D %M %L\" | grep {project}", shell=True,capture_output=True)
lines = r.stdout.decode().split("\n")
for line in lines:
    print(line)

status = np.array([x.split()[3] if len(x.split()) == 7 else "XX" for x in lines])
rtime = np.array([x.split()[6] if len(x.split()) == 7 else "0" for x in lines])
nodes = np.array([int(x.split()[4]) if len(x.split()) == 7 else 0 for x in lines])



sel = np.array(rtime != "0") * np.array(status=="PD")
hours = [int(t.split(":")[0]) + int(t.split(":")[1]) / 60 for t in rtime[sel]]
pd_time = np.sum(hours * nodes[sel])


sel = np.array(rtime != "0") * np.array(status=="R")
hours = [int(t.split(":")[0]) + int(t.split(":")[1]) / 60 for t in rtime[sel]]
run_time = np.sum(hours * nodes[sel])

sel = np.array(rtime != "0")
hours = [int(t.split(":")[0]) + int(t.split(":")[1]) / 60 for t in rtime[sel]]
tot_time = np.sum(hours * nodes[sel])
oth_time =tot_time - pd_time - run_time

# print(f"{pd_time:0.1f} SU pending for {project}")
# print(f"{run_time:0.1f} SU running for {project}")
# print(f"{oth_time:0.1f} SU other for {project}")
# print(f"{tot_time:0.1f} SU total for {project}")

print(f"{project} (assume 1:1 SU to Hours)")
print(f"Pending: \t{pd_time:>7.1f} SU")
print(f"Running: \t{run_time:>7.1f} SU")
print(f"  Other: \t{oth_time:>7.1f} SU")
print(f"--------------------------------")
print(f"  Total: \t{tot_time:>7.1f} SU")
