""""
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
rtime = np.array([x.split()[6] if len(x.split()) == 7 else "0" for x in lines])
nodes = np.array([int(x.split()[4]) if len(x.split()) == 7 else 0 for x in lines])
sel = rtime != "0"
rtime = rtime[sel]
nodes = nodes[sel]
hours = [int(t.split(":")[0]) + int(t.split(":")[1]) / 60 for t in rtime]
pd = np.sum(hours * nodes)

print(f"{pd:0.1f} SU pending for {project}")