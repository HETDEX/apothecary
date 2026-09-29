"""

usage: python avg_multifits.py <YYYYMM>

TODO: allow user to specify paths too
TODO: allow user to specify how many simultaneous procs to use

Examine averages of multi*fits for science frames
looking for issues (e.g. from --multifits_only run of reduce_shot on a month)

we are not so interested in individual shots, but in the behavior of the 78(IFU) x 4(amps) over the month

this will really need to be run in a development node ... might take longer than 2-hours ...
   may want to xlat to a simple .py file and run as a job
"""
import os.path

import numpy as np
import glob
import os.path as op
import sys
from astropy.io import fits
from tqdm  import tqdm
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file

month = int(args[0]) #202408
resume = True
SHOW_TQDM = False
MAX_WORKERS = 20 #20 works for development ... need around 8GB max to hold 300 shots of one amp in memory for
                 #averaging, etc, plus leave room for some slop
AssumedMemFootprint = 9.5 #GB (1032x1032 [pixels] x300 [shots] x8 [bytes] x2 (Python memory)

try:
    ApproxBaseRAM = psutil.virtual_memory()[0] / (1024**3) #in GB (e.g. ~32GB for vm small, 256 GB for normal on LS6)
    #need somewhere around 20GB for normal big shots (once IFU is full)
    #varies depending also on the number of exposures, but 20GB is a safe rule of thumb
    MAX_WORKERS = max(1,int(ApproxBaseRAM//AssumedMemFootprint) - 1)
    print(f"*** setting MAX_WORKERS to {MAX_WORKERS}. "
          f"BaseRAM {ApproxBaseRAM:0.1f}GB, Footprint ~ {AssumedMemFootprint}GB")
except:
    ApproxBaseRAM = -1
    MAX_WORKERS = 1


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

#todo: need to make this multithreaded ... nothing difficult here, but it will just
# take WAAAY too long otherwise ... I suppose each is independent too, so could
# do this as multiple processes, as an alternative ...


def make_avg_for_amp(amp_shot_list,amp_str):
    """
    iterate over all the shots for this one amp
    :param amp_shot_list:
    :param amp_str: (e.g. the "key" from the outer dictionary)
    :return:
    """


    if op.exists(f"{amp_str}_std.fits"):
        print(f"Skipping. Already exists: {amp_str} ...")
        return amp_str, 99

    print(f"Working on: {amp_str} ...")
    status = 0
    stop = False
    try:
        fits_list = []
        #this is the list of all multifits for that amp (e.g. the same amp for each shot in the month)
        for amp_file in tqdm(amp_shot_list,colour="green",disable=not SHOW_TQDM):

            if op.exists("stop"):
                stop = True
                break

            try:
                hdu = fits.open(amp_file)
                fits_list.append(hdu[0].data)
                hdu.close()
            except:
                print(f"Exception (inner) working on {amp_str} with file {amp_file}",traceback.format_exc())
                status = 1 #some error

        if not stop:
            keep_rt_idx = int(keep_pct * len(fits_list))

            im_array_sorted = np.sort(fits_list, axis=0)
            im_mean = np.nanmean(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
            im_median = np.nanmedian(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
            im_std = np.nanstd(im_array_sorted[0:keep_rt_idx,:,:],axis=0)

            #and/or write out as fits
            fits.PrimaryHDU(im_mean).writeto(f'{amp_str}_mean.fits', overwrite=True)
            fits.PrimaryHDU(im_median).writeto(f'{amp_str}_median.fits', overwrite=True)
            fits.PrimaryHDU(im_std).writeto(f'{amp_str}_std.fits', overwrite=True)
    except:
        print(f"Exception (outer) working on {amp_str} with file {amp_file}", traceback.format_exc())
        status = -1 #worse error

    return amp_str, status

keys = list(amp_dict.keys())
#this has one entry per unique AMP
#each iteration through the loop is independent, so this would be where to spawn other procs or make multithreaded
#main (outer) loop

with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
    futures = [pool.submit(make_avg_for_amp, amp_dict[key], key) for key in keys]

    for fut in as_completed(futures):
        if op.exists("stop"):
            break

        amp_str, amp_status = fut.result()  # re-raises any exception from the worker
        print(f"{amp_str} : status = {amp_status}")


if op.exists("stop"):
    print("stop detected. Exiting ...")
else:
    print("done")
# for key in tqdm(keys,disable=not SHOW_TQDM):
#     #fits_dict[key]     = []
#
#     if op.exists(f"{key}_std.fits"):
#         print(f"Skipping. Already exists: {key} ...")
#         continue
#
#
#


    # print(f"Working on: {key} ...")
    # fits_list = []
    # #this is the list of all multifits for that amp (e.g. the same amp for each shot in the month)
    # for amp in tqdm(amp_dict[key],colour="green",disable=not SHOW_TQDM):
    #     hdu = fits.open(amp)
    #     #fits_dict[key].append(hdu[0].data)
    #     fits_list.append(hdu[0].data)
    #
    # #keep_rt_idx = int(keep_pct * len(fits_dict[key]))
    # keep_rt_idx = int(keep_pct * len(fits_list))
    #
    # #im_array_sorted = np.sort(fits_dict[key],axis=0)
    # im_array_sorted = np.sort(fits_list, axis=0)
    # im_mean = np.nanmean(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    # im_median = np.nanmedian(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    # im_std = np.nanstd(im_array_sorted[0:keep_rt_idx,:,:],axis=0)
    # # mean_dict[key] = im_mean
    # # median_dict[key] = im_median
    # # std_dict[key] = im_std
    #
    # #and/or write out as fits
    # fits.PrimaryHDU(im_mean).writeto(f'{key}_mean.fits', overwrite=True)
    # fits.PrimaryHDU(im_median).writeto(f'{key}_median.fits', overwrite=True)
    # fits.PrimaryHDU(im_std).writeto(f'{key}_std.fits', overwrite=True)

#todo: or in another python code,
# basic analysis of each averaged mutli*fits
#   highlight (color?) pixels that are high and low (in mean or median averages) where
#          high/low is a sigma? or some fixed exteme value?
#  highlight (color?) pixels that are high in std fits

# return a report that identifies amps that have many such pixels for manual review
# amps that just have a few are maybe ignorable? or perhaps have the report list the amp and the
#    total numbers of "suspect" pixels and sort by those totals so we check on the worst ones first?