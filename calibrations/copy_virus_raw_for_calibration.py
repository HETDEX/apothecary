"""

for a month, copy the raw files needed for virus calibrations

these either come from /work/03946/hetdex/maverick/YYYYMM
or they come from /corral/utexas/Hobby-Eberly-Telesco/het_raw/YYYYMM

rather than copy the entire directory (maverick) or the whole tar (corral), just copy to
your own het_raw the bits you need

These would be mostly the same as in reduce_shot.py, the gc1, gc2, and the virus0000XXX shot tars

input is a YYYYMM

output is a set of folders in the current working directory, one for each day of the YYYYMM

"""


import sys
import os
import glob
from pathlib import Path
import numpy as np
from dataclasses import dataclass
import traceback
import tarfile as tar
from tqdm import tqdm



args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]


if len(args) > 0:
    yyyymm = int(args[0])
else:
    print("Must specify yyyymm")
    exit(-1)

if not 201601 <= yyyymm <=209912:
    print(f"Invalid yyyymm {args[0]}")
    exit(-1)

#HETRaw_archives = ["/work/03946/hetdex/maverick","/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"]
corral_basedir = "/corral/utexas/Hobby-Eberly-Telesco/het_raw"
maverick_basedir = "/work/03946/hetdex/maverick"  #note calling it maverick to not confuse with my own "work"
# print("TESTING!!!! Alternate corral_basedir, maverick_basedir")
# corral_basedir = "/scratch/03261/polonius/het_raw_test_src_corral"
# maverick_basedir = "/scratch/03261/polonius/het_raw_test_src_maverick"
het_raw_path = None



##################################
# pseudo MAIN
##################################


#build list of each and compbne to be unique with a prefernce for the already "exploded" directories on maverick
fns_maverick = glob.glob(os.path.join(maverick_basedir,f"{yyyymm}*"))
fns_corral = glob.glob(os.path.join(corral_basedir,f"{yyyymm}*.tar"))

dates_maverick = [int(os.path.basename(fn).split(".")[0]) for fn in fns_maverick]
dates_corral = [int(os.path.basename(fn).split(".")[0]) for fn in fns_corral] #strip off the .tar, just keep YYYYMMDD

#combine for unique list
#where there are duplicates, keep maverick
#basically iterate over maverick (which could be empty) if it has a match in corral, delete the corral match from its list
for date in dates_maverick:
    try:
        j = dates_corral.index(date)
        del fns_corral[j]
        del dates_corral[j]
    except:
        pass #not in corral


#now we have two lists without duplicates


#need subdirs gc1, gc2, virus  and nothing else (no acm, dimm, hpf, or lsr2)
#each subdir has one or more *.tar files

if het_raw_path is None: #if not configured, you are copying into the current directory
    het_raw_path = "./"

#could do maverick then corral OR go by date order, flipping back and forth?
if len(fns_maverick) > 0:
    print(f"Copy from {maverick_basedir} ... ")
    for fn in tqdm(fns_maverick): #top level is a folder for the date
        #do the copy
        try:
            yyyymmdd = os.path.basename(fn) #just a folder
            dest = os.path.join(het_raw_path, yyyymmdd)
            if not os.path.exists(dest):
                Path(dest).mkdir(parents=True, exist_ok=True)

            #need subdirs gc1, gc2, virus  and nothing else? (no acm, dimm, hpf, or lsr2
            #do not copy if already there
            src = os.path.join(fn,"gc1")
            if not os.path.exists(os.path.join(dest,"gc1/gc1.tar")):
                cmd = f"cp -r {src} {dest}"
                os.system(cmd)
            #else:
            #    print("skipping gc1")

            src = os.path.join(fn,"gc2")
            if not os.path.exists(os.path.join(dest, "gc2/gc2.tar")):
                cmd = f"cp -r {src} {dest}"
                os.system(cmd)
            #else:
            #    print("skipping gc2")

            src = os.path.join(fn,"virus")
            subfiles_dest = sorted(glob.glob(f"{dest}/virus/virus*.tar"))

            if len(subfiles_dest) == 0:
                cmd = f"cp -r {src} {dest}"
                os.system(cmd)
            else:
                #some might be here
                subfiles_src = sorted(glob.glob(f"{src}/virus*.tar"))
                fns_src = [os.path.basename(s) for s in subfiles_src]
                fns_dest = [os.path.basename(s) for s in subfiles_dest]
                sel_src = np.array([s not in fns_dest for s in fns_src])
                subfiles_src = np.array(subfiles_src)[sel_src]

                #print(f"Copying limited number of files: {subfiles_src}")

                #dest += "/virus"
                for fn in subfiles_src:
                    cmd = f"cp -r {src} {dest}"
                    os.system(cmd)

        except:
            print(f"[{yyyymmdd}] Failed to copy date dir.", traceback.format_exc())


if len(fns_corral) > 0:
    print(f"Copy/Untar from {corral_basedir} ... ")
    for fn in tqdm(fns_corral): #top level is a tar for the date
        # this is the date.tar, then need to extract parts of it
        try:
            yyyymmdd = os.path.basename(fn).split(".")[0] #strip off the ".tar

            # if not os.path.exists(os.path.join(het_raw_path,yyyymmdd)):
            #     Path(os.path.join(het_raw_path,yyyymmdd)).mkdir(parents=True, exist_ok=True)

            file_to_extract = f"{yyyymmdd}/gc1"
            cmd = f"tar -xvf {fn} --skip-old-files -C {het_raw_path} {file_to_extract}"
            os.system(cmd)

            file_to_extract = f"{yyyymmdd}/gc2"
            cmd = f"tar -xvf {fn} --skip-old-files -C {het_raw_path} {file_to_extract}"
            os.system(cmd)

            file_to_extract = f"{yyyymmdd}/virus"
            cmd = f"tar -xvf {fn} --skip-old-files -C {het_raw_path} {file_to_extract}"
            os.system(cmd)

        except:
            print(f"[{yyyymmdd}] Failed to copy/extract date tar file.", traceback.format_exc())


