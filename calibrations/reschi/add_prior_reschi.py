"""
related to and to be run after make_reshi.py

This is the part that takes the prior reschi and adds it in to the current (newly generated) resXXXAA.fits files

This is run under the parallel reduction directory, normally under the directory for a specific month, and under
   that directory's "reductions" folder (e.g. <somepath>/parallel/m202408/reductions/)
   The directory above should have all the sci<YYYYMMDDvSSS> folders


"""
import glob
import os
import numpy as np
from tqdm import tqdm
from astropy.io import fits
import traceback

#BASEDIR_reschi = "/scratch/projects/hetdex/lib_calib/reschi" #do not end with /
#to use local reschi, copy the per specID contents to reschi/yyyymm/ locally and set BASEDIR_reschi to that path
BASEDIR_reschi = "/scratch/03261/polonius/parallel/reschi" #do not end with /


print(f"Using reschi under: {BASEDIR_reschi}")

all_dateshot_paths = sorted(glob.glob("20*/virus/virus*"))
all_dates = [x.split("/")[0] for x in all_dateshot_paths]
all_shotnums = [x.split("/")[-1][-3:] for x in all_dateshot_paths]


# !!!  assume the immediately preceeding month is the prior
#      if this is not so, edit here !!!
current_year = int(all_dates[0][:4]) #string YYYY
current_month = int(all_dates[0][4:6]) #string MM
current_yyyymm = all_dates[0][:6]
if current_month == 1:
    prior_month = 12
    prior_year = current_year - 1
else:
    prior_year = current_year
    prior_month = current_month -1

prior_yyyymm = f"{prior_year}{str(prior_month).zfill(2)}"

#I do NOT want these sorted ... the two arrays need to be in the same order
all_prior_res_files = glob.glob(f"{BASEDIR_reschi}/{prior_yyyymm}/res*.fits")
all_prior_res_files_basenames = [os.path.basename(x) for x in all_prior_res_files]

#get all the currently generate fits files
all_new_res_files = glob.glob(f"spec*/res*.fits")
all_new_res_files_basenames = [os.path.basename(x) for x in all_new_res_files]


#match them up, and add
ct_already_updated = 0
ct_newly_updated = 0
ct_failed = 0
list_failed = [] #basically specids in the current month that are not in the prior month

for new_fn in tqdm(all_new_res_files_basenames):
    new_ix = None
    prior_ix = None
    try:
        new_ix = all_new_res_files_basenames.index(new_fn)
        if os.path.exists(f"{all_new_res_files[new_ix]}.done"):
            #this one was already done
            #print(f"{new_fn} already updated with prior reschi")
            ct_already_updated += 1
            continue

        #find the matching index of the prior files
        prior_ix = all_prior_res_files_basenames.index(new_fn)

        #normally, these would be at the same index, but I suppose that is not strictly true

        prior_hdu = fits.open(all_prior_res_files[prior_ix])
        new_hdu = fits.open(all_new_res_files[new_ix],mode='update')

        new_hdu[0].data += prior_hdu[0].data
        new_hdu.flush()

        new_hdu.writeto(all_new_res_files[new_ix], overwrite=True)

        with open(f"{all_new_res_files[new_ix]}.done", "w") as f:
            f.write(f"Added in {all_prior_res_files[prior_ix]}")
            ct_newly_updated += 1

        new_hdu.close()
        prior_hdu.close()

    except:
        #print(f"Unable to add prior reschi to {new_fn}: {traceback.format_exc()}")
        if prior_ix is None:
            #print(f"Unable to add prior reschi to {new_fn}, no matching file")
            pass #will show up in later count
        else:
            print(f"Unable to add prior reschi to {new_fn}: {traceback.format_exc()}")
        ct_failed += 1
        list_failed.append(new_fn)



#we only care about mismatches for this purpose
#files in prior month that are not in the current month
list_prior_month_no_match = np.setdiff1d(all_new_res_files_basenames,all_prior_res_files_basenames)
#files in the current month that are not in the prior month
list_current_month_no_match = np.setdiff1d(all_prior_res_files_basenames,all_new_res_files_basenames)

print(f"New updates : {ct_newly_updated}")
print(f"Already done: {ct_already_updated}")
print(f"Failed      : {ct_failed}")
for fn in sorted(list_failed):
    if fn in list_prior_month_no_match:
        print(f"  {fn} not in prior month")
    else:
        print(f"  {fn} failed - unknown")

print(f"\nCurrent specIDs that are not in prior month: \n{list_prior_month_no_match}")
print(f"\nPrior specIDs that are not in current month: \n{list_current_month_no_match}")
