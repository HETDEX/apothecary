"""
!!!! DO NOT USE THIS ONE!!!! UNLESS YOU KNOW WHAT YOU ARE DOING !!!!

Generally use get_rescor_by_shot_nomask

"""



import numpy as np
import sys
import os
import h5py

import traceback

from astropy.table import Table

from hetdex_api.shot import open_shot_file, get_fibers_table

from hetdex_lim.ifu_spectra import IFU_Spectra
from hetdex_lim.config import LIM_Config

HDRVersion = "hdr5"

if __name__ == "__main__":
    no_overwrite = True #do not re-run if all files for shot have been done
    version = sys.argv[1]
    shotid = int(sys.argv[2])
    month = str(shotid)[:-5]
#    save_filename_dir = f"/scratch/05865/maja_n/ffskysub/rescor/mask/{month}/individual/"
    save_filename_dir = f"/scratch/03261/polonius/rescor/code/rescor/ffskysub/rescor/mask/{month}/individual/"
    if not os.path.exists(save_filename_dir):
        try:
            print(f"mkdir {save_filename_dir}")
            os.makedirs(save_filename_dir)
        except Exception as e:  #with multiples running this could get created by another
            if not os.path.exists(save_filename_dir):
                print(f"Problem. Exception attempting to create output dir {save_filename_dir}")
                exit(-1)


    save_filename_pattern = os.path.join(save_filename_dir, "{}_{}_{}.h5")

    lim_config = LIM_Config(version)
    def_wave = np.arange(3470, 5542, 2)

    ifuspec = IFU_Spectra(shotid,
                        save_filename = os.path.join(lim_config.shot_files_dir, f"{shotid}.h5"),
                            do_boxcar = lim_config.do_boxcar,
                            do_mask_crazy_pixels = lim_config.do_mask_crazy_pixels,
                            do_mask_continuum = lim_config.do_mask_continuum,
                            do_mask_LAEs = lim_config.do_mask_LAEs,
                            LAE_catalog = lim_config.lae_catalog_filename,
                            do_split_fibers = lim_config.do_split_fibers,
                            fiber_mask_filename = lim_config.fiber_mask_filename,
                            sncut = lim_config.sncut)

    ifuspec.load_shot()
    ifuspec.mask_badamps()
    ifuspec.mask_bad_fibers()
    ifuspec.mask_continuum_emission()

    try:
        fileh = open_shot_file(shotid,survey=HDRVersion)
        ffskymod_1 = fileh.root.FullSkyModel.exp01.read()
        ffskymod_2 = fileh.root.FullSkyModel.exp02.read()
        ffskymod_3 = fileh.root.FullSkyModel.exp03.read()
        fileh.close()
    except Exception as e:
        print(f"Exception loading shot file for {shotid}. Must abort.")
        print(f"{shotid} Exception: {e}\n\n{traceback.format_exc()} ")
        exit(-1)

    ffskymod = {}
    for exp_num, model in enumerate([ffskymod_1, ffskymod_2, ffskymod_3]):
        wave, spec = np.transpose(model) # spec is in counts
        ffskymod[exp_num + 1] = np.interp(def_wave, wave, spec)
    del ffskymod_1, ffskymod_2, ffskymod_3

    multiframes = np.unique(ifuspec.multiframe)
    exps = np.unique(ifuspec.expnum)


    should_continue = False
    if no_overwrite:
        for expnum in exps:
            for multiframe in multiframes:
                save_filename = save_filename_pattern.format(shotid, multiframe, expnum)
                if not os.path.isfile(save_filename):
                    should_continue = True

    if not should_continue:
        print(f"All multiframes and exposures already done for this shot {shotid}")
    else:
        fibers = get_fibers_table(shotid,survey=HDRVersion)
        gain = fibers['calfib'] / fibers['calfib_counts']
        del fibers

        for expnum in exps:
            for multiframe in multiframes:
                here = (ifuspec.expnum==expnum) * (ifuspec.multiframe==multiframe)

                res_here = ifuspec.ffskysub[here] / gain[here] / ffskymod[expnum][None,:]

                save_filename = save_filename_pattern.format(shotid, multiframe, expnum)

                with h5py.File(save_filename, "w") as ff:
                    ff['rescor'] = res_here
