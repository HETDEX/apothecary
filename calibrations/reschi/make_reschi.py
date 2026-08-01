"""
make the reschi (res*.fits and rres*.fits) for a month


This is run under the parallel reduction directory, normally under the directory for a specific month, and under
   that directory's "reductions" folder (e.g. <somepath>/parallel/m202408/reductions/)
   The directory above should have all the sci<YYYYMMDDvSSS> folders

This is part 1 of 3. This will produce rgetres.sh and rgetres.run. This is very fast (just seconds).
As part 2 of 3, run either ... rgetrsh.sh in an idev session (vm-small is fine and it can take more than 1 hour)
   or build a slurm and run rgetres.run in that slurm.
After this runs, part 3 of 3 is to run add_prior_reschi.py



this is a wrapper for Karl's rgetres0, rgetres1 which call in to ~/gebhardt/bin/


Notes: from email 20260627
res*fits : average residual after sky subtraction (112x1036 arrays). This is designed to provide an additional
sky subtraction in vred. It is important to get this correct. This is only used in vred.

rres*fits : stands for ratio of the residuals, which is the error at that element. The idea here was to use as an
additional uncertainty when measuring s/n. This was too hard to get correct, so I dropped it. Currently it is only used
 as a flag in fitradecsp, and removes elements that have the value<-0.05

chi*fits: not really informative, but it looks at the distribution of residuals at that element and measure a chi^2

Well, turns out that none of those arrays are being used, except for in extreme cases. Looks like I turned all this off
 some time ago (years ago) since the idea was to get the sky subtraction correct from the start.

The only thing that is being used is if rres is below -0.5, then that element is flagged to not be used.

It is coming back to me since the proper thing to do is to fix it. I just could never trust the analysis of the
residual corrections for detections.

So, fitradecsp does not need any of them, except for flags. However ...
I had put the residual correction in the main reductions of vred. And that is pointing to the old setup.
For that correction, it only uses the res*.fits

Ok, so this changes what we do. Let's only track res by month, and use rres for flags.
So then each month in reschi would have a full suite of res*fits, and then only a handful of rres*fits that have flags.


Maybe in a future date, we can address ability to use rres as intended. We might be able to explore with some of the
high sky frames.

This explains why it was not doing what I thought it should be doing.


basic logic

for a given month,
1. run reduce_shot --multifits_only YYYYMM     for that month
2. in the output reductions folder, do some prep work
   * get a list of all unique SPECIDs
     + this becomes the calls to rgetres1
   * get a list of all shots (date and shot number, e.g. datevshot)
     + turn this into a file called "mlist" which is a series of lines as:  YYYYMMDD SSS
   * get a list of all res*.fits for the previous month *** (assumes this is the active reschi for the reduce_shot call)
     + this will be used at the end to add into the initial resulting res*.fits files
3.

"""

import glob
import os
import numpy as np

BASEDIR_reschi = "/scratch/projects/hetdex/lib_calib/reschi" #do not end with /

#Assume reduce_shot --multifits_only YYYYMM has already completed and you are in the reductions directory
#e.g. /scratch/03261/polonius/parallel/m202408/reductions/20240801/virus/virus0000002/exp01/virus
#find all the unique SPECIDs for exp01 only
all_fits_fns = glob.glob("20*/virus/virus*/exp01/virus/multi_*.fits")
all_fits_fns = [os.path.basename(x) for x in all_fits_fns]
unique_fits_fns = np.unique(all_fits_fns)
unique_specids = sorted(np.unique([x.split("_")[1] for x in unique_fits_fns]))


all_dateshot_paths = sorted(glob.glob("20*/virus/virus*"))
all_dates = [x.split("/")[0] for x in all_dateshot_paths]
all_shotnums = [x.split("/")[-1][-3:] for x in all_dateshot_paths]

current_yyyymm = all_dates[0][:6]

# moved to add_prior_reschi.py
# # !!!  assume the immediately preceeding month is the prior
# #      if this is not so, edit here !!!
# current_year = int(all_dates[0][:4]) #string YYYY
# current_month = int(all_dates[0][4:6]) #string MM
# current_yyyymm = all_dates[0][:6]
# if current_month == 1:
#     prior_month = 12
#     prior_year = current_year - 1
# else:
#     prior_year = current_year
#     prior_month = current_month -1
#
# prior_yyyymm = f"{prior_year}{str(prior_month).zfill(2)}"
#
# all_prior_res_files = sorted(glob.glob(f"{BASEDIR_reschi}/{prior_yyyymm}/res*.fits"))
# all_prior_res_files = [os.path.basename(x) for x in all_prior_res_files]

######################################
# rgetres1 prep
######################################
with open("mlist","w") as f:
    for date, shot in zip(all_dates, all_shotnums):
        f.write(f"{date} {shot}\n")

#this would be for a slurm
with open("rgetres.run","w") as f:
    for specid in unique_specids:
        f.write(f"./rgetres1 {specid} {current_yyyymm}\n")

#this would be to run locally in an idev sesssion
#normally there should be 78 calls, lets make that 10 calls with 8 per line (last line is short)
with open("rgetres.sh","w") as f:
    for i, specid in enumerate(unique_specids):

        if i % 8 == 0: #this is a first line
            if i != 0:
                f.write(" ) & \n") #close off the last line

            f.write(f"( ./rgetres1 {specid} {current_yyyymm} ")
        else:
            f.write(f" && ./rgetres1 {specid} {current_yyyymm} ")

    if i % 8:
        f.write(" ) & \n")  # close off the last line
    f.write("wait\n")

#now run the reschi generation #maybe a limited slurm?

print("Now open an idev and run rgetres.sh OR run as a SLURM.")
print("After that is done, do the addition of the previous month.")


####################################################
#Now need to do the prior month adding ... this should be in another Python file
####################################################
