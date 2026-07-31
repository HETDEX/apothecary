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


corral_basedir = "/corral/utexas/Hobby-Eberly-Telesco/het_raw"
maverick_basedir = "/work/03946/hetdex/maverick"  #note calling it maverick to not confuse with my own "work"
het_raw_path = None

#build list of each and compbne to be unique with a prefernce for the already "exploded" directories on maverick
fns_maverick = glob.glob(os.path.join(maverick_basedir,f"{yyyymm}*"))
fns_corral = glob.glob(os.path.join(corral_basedir,f"{yyyymm}*.tar"))

dates_maverick = [int(os.path.basename(fn)) for fn in fns_maverick]
dates_corral = [int(os.path.basename(fn)[0:8]) for fn in fns_corral] #strip off the .tar, just keep YYYYMMDD

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

#could do maverick then corral OR go by date order, flipping back and forth?
for fn in fns_maverick:
    #do the copy
    pass

for fn in fns_corral:
    #copy from tar (see reduce_shot.py)
    pass

