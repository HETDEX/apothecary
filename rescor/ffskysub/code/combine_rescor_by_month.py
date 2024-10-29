# run example (January 2022)
# python combine_rescor_by_month.py 202201 1 0

import numpy as np
import sys
import os
import glob
import h5py
import copy
import traceback
from multiprocessing import Pool

from astropy.stats import biweight_location
from astropy.io import fits
from astropy.table import Table

from hetdex_api.survey import Survey
from hetdex_api.config import HDRconfig

def combine_multiframe(multiframe):

    """
    This function takes the multiframe 
    and uses global variables like files 
    and save_filename_pattern.

    It opens all rescor files in 'files',
    masks zeros -> NaN
    calculates the biweight location (as a function of fiber and wavelength)
    and saves this biweight rescor to an h5 file.
    """

    here = multiframes == multiframe
    all_rescors = []
    for filename in files[here]:
        with h5py.File(filename, "r") as ff:
            if ff['rescor'].size == 0:
                continue
            all_rescors.append(ff['rescor'][:])
    all_rescors = np.array(all_rescors)
    N_files = len(all_rescors)
    print(f"{all_rescors[0].shape=}")
    all_rescors[all_rescors == 0] = np.nan
    av_rescor = biweight_location(all_rescors, axis=0, ignore_nan=True)

    with h5py.File(save_filename_pattern.format(multiframe), "w") as ff:
        ff['rescor'] = av_rescor
        ff['rescor'].attrs['N_files'] = N_files

    print(f"Done with {multiframe}; {save_filename_pattern.format(multiframe)}")
    return


if __name__ == '__main__':  

    # input #1 can be a month (recommended) or a year 
    # input #2 is 1 if it only masks bad amps (recommended)
    # input #3 is 0 if you do not split by position angle (recommended)
    month = str(sys.argv[1])
    nomask = int(sys.argv[2]) == 1
    split_by_pa = int(sys.argv[3]) == 1

    all_in_one_dir = True #all the shot.h5 are in a single directory
    
    # check if 'month' variable is a month or a year.
    if len(month) == 4:
        # then this is actually a year and we should combine the files of the whole year. 
        do_by_month = False
    elif len(month) == 6:
        # then this is a month and we combine all files in one month directory.
        do_by_month = True
    do_by_year = ~do_by_month

    # if we split by PA, check the PA of all shots
    if split_by_pa:
        survey = Survey("hdr4").return_astropy_table()
        pa_group_1 = {survey['shotid'][i] : survey['pa'][i] < 180 for i in range(len(survey))}
        pa_group_2 = {survey['shotid'][i] : (survey['pa'][i] >= 180) for i in range(len(survey))}

    # get the file directory of applicable rescor files (depending on mask/nomask)
    #operating_basedir = "/scratch/03261/polonius/rescor/code/rescor/ffskysub/code"
    operating_basedir = "./"
    #read_filename_dir = "/scratch/05865/maja_n/ffskysub/rescor/"
    #read_filename_dir = "/scratch/03261/polonius/red1/rescor/ffskysub/rescor/"
    #read_filename_dir = "/scratch/03261/polonius/rescor/code/rescor/ffskysub/improved_spectra/pytables/"
    read_filename_dir = os.path.join(operating_basedir,"./output/ffskysub/rescor/")
    if nomask:
        read_filename_dir = os.path.join(read_filename_dir, "nomask")
    else:
        read_filename_dir = os.path.join(read_filename_dir, "mask")
    
    # get the pattern of rescor files depending on whether to average by year or by month
    #if all_in_one_dir:
    #    print("all in one")
    #    read_filename_pattern = os.path.join(read_filename_dir, f"{month}*.h5")
    if do_by_month:
        print("By month")
        read_filename_pattern = os.path.join(read_filename_dir, f"{month}/individual/2*.h5")
    elif do_by_year:
        print("By year")
        read_filename_pattern = os.path.join(read_filename_dir, f"{month}??/individual/2*.h5")


    # define the directory for the output files
    #save_filedir_0 = "/".join(read_filename_dir.split("/")[3:])

    save_filedir_0 = copy.copy(read_filename_dir)
    save_filedir = os.path.join(save_filedir_0, month, "combined")
    if split_by_pa:
        save_filedirs = [os.path.join(save_filedir, f"pa_{idx}") for idx in ['E', 'W']]
    else:
        save_filedirs = [save_filedir]

    # define the pattern for the output files and make sure that the output directories exist
    save_filename_patterns = [os.path.join(save_filedir, "combined_{}.h5") for save_filedir in save_filedirs]
    for save_filedir in save_filedirs:
        print(f"checking {save_filedir}")
        if not os.path.exists(save_filedir):
            print(f"mkdir {save_filedir}")
            os.makedirs(save_filedir)
        else:
            print("exists")

    # get all files
    all_files = np.array(sorted(glob.glob(read_filename_pattern)))

    #if there are not enough, get the previous and next months? and append to all_files
    #what is enough? say. minimum of 30x3 exposures (so 90) per multiframe, meaning 30 shots
    #so ... on a per multiframe basis, if there are < 90, grab those from the previous and next month, if any
    if do_by_month:
        add_flanking_months = False
        try:
            #strip off just th e IFU and AMP bit,
            mf, ct = np.unique([s[-19:-5] for s in all_files],return_counts=True)
            #todo: *could* do this by multiframe ... so some might pull from last and next month and some might not

            sel = ct < 90
            add_flanking_mfs = mf[sel]
            if np.count_nonzero(sel):
                add_flanking_mfs = mf[sel]
                add_flanking_months = True

        except Exception as E:
            print(f"Exception checking all_files: {e}\n\n{traceback.format_exc()}")

        if add_flanking_months:
            print(f"Adding flanking months in for {len(add_flanking_mfs)} multiframes")
            last_month = np.datetime64(f"{month[0:4]}-{month[4:]}") - np.timedelta64(1,"M")
            next_month = np.datetime64(f"{month[0:4]}-{month[4:]}") + np.timedelta64(1,"M")

            last_month = str(last_month).replace('-', '')
            next_month = str(next_month).replace('-', '')

            for mf in add_flanking_mfs:
                #last
                pattern = os.path.join(read_filename_dir, f"{last_month}/individual/2*_multi_{mf}_?.h5")
                flanking_files = np.array(sorted(glob.glob(pattern)))
                all_files = np.concatenate((all_files,flanking_files))

                #next
                pattern = os.path.join(read_filename_dir, f"{next_month}/individual/2*_multi_{mf}_?.h5")
                flanking_files = np.array(sorted(glob.glob(pattern)))
                all_files = np.concatenate((all_files, flanking_files))


    # if split by PA, split the files by PA
    if split_by_pa:
        shotids = [int(ff.split("/")[-1].split("_")[0]) for ff in all_files]
        in_group_1 = [pa_group_1[shotid] for shotid in shotids]
        in_group_2 = [pa_group_2[shotid] for shotid in shotids]
        len_all_files = len(all_files)
        all_files = [all_files[in_group_1], all_files[in_group_2]]
        assert len_all_files == np.sum([len(all_files[i]) for i in range(len(all_files))])
    else:
        all_files = [all_files]

    # get badamp dictionary
    # this is redundant if you ran get_rescor_by_shot_nomask.py 
    # in the version that included the bad amp masking.
    # if there might still be old files from bad amps, 
    # this is necessary to exclude them.
    if True:
        print("Getting badamp dict")
        config = HDRconfig()
        badamp_tab  = fits.open(config.badamp)
        badamp_tab = Table(badamp_tab[1].data)
        # columns: 'shotid', 'multiframe', 'flag'
        badamp_dict = {(badamp_tab['shotid'][i], badamp_tab['multiframe'][i]) : badamp_tab['flag'][i] for i in range(len(badamp_tab))}
        print("Got badamp dict.")
    else:
        print("Skipping badamp dict")
        badamp_dict = {}
    
    # if split by PA, loop through the PAs. 
    # otherwise, all_files and save_filename_patterns have length 1
    ## DD appears that all_files is expected to contain names like: 20221023022_multi_414_038_035_RL_3.h5
    ## as might be under the "YYYYMM/individual/YYYYMMDDSSS_multi*h5"
    for save_filename_pattern, files in zip(save_filename_patterns, all_files):
        multiframes = ["_".join(x.split("/")[-1].split("_")[1:-1]) for x in files]
        shotids = [int(x.split("/")[-1].split("_")[0]) for x in files]
        files = list(files)

        # remove bad amps from the file list 
        i = 0
        while i < len(shotids):
            try:
                if badamp_dict[(shotids[i], multiframes[i])] == 0:
                    print(f"Removing bad amp {shotids[i]} {multiframes[i]}")
                    del files[i], shotids[i], multiframes[i]
                else:
                    i += 1
            except:
                i += 1

        # get unique multiframes
        un_multiframes = list(np.unique(multiframes))

        # if this is set to True, it will not re-calculate and overwrite 
        # combined rescor files that already exist.
        if False:
            i = 0
            while i < len(un_multiframes):
                if os.path.exists(save_filename_pattern.format(un_multiframes[i])):
                    del un_multiframes[i]
                else:
                    i += 1

        multiframes = np.array(multiframes)
        files = np.array(files)

        # combine the multiframes for each multiframe
        # where it is not a bad amp.
        if len(un_multiframes) > 0:
            with Pool(8) as pool:
                results = pool.map(combine_multiframe, un_multiframes)