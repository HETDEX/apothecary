"""
related to and to be run after make_reshi.py

This is the part that takes the prior reschi and adds it in to the current (newly generated) resXXXAA.fits files


"""
import glob
import os
import numpy as np
from tqdm import tqdm
from astropy.io import fits
import traceback

BASEDIR_reschi = "/scratch/projects/hetdex/lib_calib/reschi" #do not end with /


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
all_new_res_files_basenames = [os.path.basename(x) for x in all_prior_res_files]


#match them up, and add
ct_already_updated = 0
ct_newly_updated = 0
ct_failed = 0
for new_fn in tqdm(all_new_res_files_basenames):
    try:
        new_ix = all_new_res_files_basenames.index(new_fn)
        if os.path.exists(f"{all_new_res_files[new_ix]}.done"):
            #this one was already done
            print(f"{new_fn} already updated with prior reschi")
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
        print(f"Unable to add prior reschi to {new_fn}: {traceback.format_exc()}")
        ct_failed += 1

print(f"New updates : {ct_newly_updated}")
print(f"Already done: {ct_already_updated}")
print(f"Failed      : {ct_failed}")