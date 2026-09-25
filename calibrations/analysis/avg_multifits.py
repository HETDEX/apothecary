"""
Examine averages of multi*fits for science frames
looking for issues (e.g. from --multifits_only run of reduce_shot on a month)

we are not so interested in individual shots, but in the behavior of the 78(IFU) x 4(amps) over the month

this will really need to be run in a development node ... might take longer than 2-hours ...
   may want to xlat to a simple .py file and run as a job
"""


import numpy as np
import glob
import os.path as op
from astropy.io import fits
from tqdm  import tqdm

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


month = 202408
resume = True



reductions_basepath = f"/scratch/03261/polonius/parallel/m{month}" #basically gets the month
#from there the paths are YYYYMMDD/virus/virus0000<shot>/exp<##>/virus/multi*
all_multis = np.array(glob.glob(f"{reductions_basepath}/reductions/*/virus/virus0000*/exp*/virus/multi*.fits"))
all_amps = np.array([op.basename(x).split(".")[0] for x in all_multis]) #need this as an array so can select on it
unique_amps = np.unique(all_amps)


#organize in dictionary so can deal with each amp

print(f"Building averages of multi*fits for {month} from {reductions_basepath} ...")
amp_dict = {}
for amp in unique_amps:
    sel = all_amps == amp
    amp_dict[amp] = all_multis[sel]


#fits_dict = {}
#mean_dict = {}
#median_dict = {}
#std_dict = {}

keep_pct = 0.95

keys = list(amp_dict.keys())
for key in tqdm(keys):
    #fits_dict[key]     = []
    fits_list = []
    if op.exists(f"{key}_std.fits"):
        print(f"Skipping. Already exists: {key} ...")
        continue

    print(f"Working on: {key} ...")
    for amp in tqdm(amp_dict[key],colour="green"):
        hdu = fits.open(amp)
        #fits_dict[key].append(hdu[0].data)
        fits_list.append(hdu[0].data)

    #keep_rt_idx = int(keep_pct * len(fits_dict[key]))
    keep_rt_idx = int(keep_pct * len(fits_list))

    #im_array_sorted = np.sort(fits_dict[key],axis=0)
    im_array_sorted = np.sort(fits_list, axis=0)
    im_mean = np.nanmean(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    im_median = np.nanmedian(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    im_std = np.nanstd(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    # mean_dict[key] = im_mean
    # median_dict[key] = im_median
    # std_dict[key] = im_std

    #and/or write out as fits
    fits.PrimaryHDU(im_mean).writeto(f'{key}_mean.fits', overwrite=True)
    fits.PrimaryHDU(im_median).writeto(f'{key}_median.fits', overwrite=True)
    fits.PrimaryHDU(im_std).writeto(f'{key}_std.fits', overwrite=True)

#todo: or in another python code,
# basic analysis of each averaged mutli*fits
#   highlight (color?) pixels that are high and low (in mean or median averages) where
#          high/low is a sigma? or some fixed exteme value?
#  highlight (color?) pixels that are high in std fits

# return a report that identifies amps that have many such pixels for manual review
# amps that just have a few are maybe ignorable? or perhaps have the report list the amp and the
#    total numbers of "suspect" pixels and sort by those totals so we check on the worst ones first?