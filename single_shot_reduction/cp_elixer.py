"""
  helper to copy elixer catalog h5 and all_pngs to staging area on /scratch/projects
  creates a batch file to be sourced as the hetdex user

  this py file needs to be run in the elixer directory
  and the resulting bash file in the target direcetory

  this takes a few minutes per elixer (the cost is in the all_pngs files) ... is it worth zipping? we are just copying
    from /scratch to /scratch though ... same filesytem

"""


import glob
import os
#from pathlib import Path

h5_fns = sorted(glob.glob("sci*/elixer/out/elixer_merged_cat.h5"))
basedir = os.getcwd()

with open("cp_elixer","w") as f:
    for fn in h5_fns:
        topdir = fn.split("/")[0]
        cath5 = os.path.join(basedir,fn)
        allpngs = os.path.join(basedir,topdir,"elixer/out/all_pngs")

        cmdstr = f"mkdir -p {topdir} ; cd {topdir} ; cp {cath5} . ; cp -r {allpngs} . ; cd .. \n"

        f.write(f"echo {topdir} ... \n")
        f.write(cmdstr)


