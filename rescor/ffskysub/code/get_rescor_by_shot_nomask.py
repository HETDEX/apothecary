# run example:
# python get_rescor_by_shot_nomask.py 20220103007

import numpy as np
import sys
import os
import h5py

import traceback

from astropy.table import Table
from astropy.io import fits

from hetdex_api.shot import open_shot_file, get_fibers_table
from hetdex_api.config import HDRconfig


HDRVersion = "hdr5"

if __name__=="__main__":

    no_overwrite = True  # do not re-run if all files for shot have been done

    # the input is the shotid
    shotid = int(sys.argv[1])
    month = str(shotid)[:-5]

    #working_dir = "/scratch/03261/polonius/rescor/code/rescor/ffskysub/code"
    working_dir = "./"
    # make sure that the directory to save files to exists
#    save_filename_dir = f"maja_n/ffskysub/rescor/nomask/{month}/individual/"
    save_filename_dir = os.path.join(working_dir,f"output/ffskysub/rescor/nomask/{month}/individual/")
    if not os.path.exists(save_filename_dir):
        print(f"mkdir {save_filename_dir}")
        os.makedirs(save_filename_dir)
    save_filename_pattern = os.path.join(save_filename_dir, "{}_{}_{}.h5")

    running_status_file = os.path.join(save_filename_dir,f"rescor_by_shot_nomask_{shotid}.running")
    done_status_file = os.path.join(save_filename_dir, f"rescor_by_shot_nomask_{shotid}.done")

    try:
        if no_overwrite and os.path.isfile(done_status_file):
            print(f"{shotid} All multiframes and exposures already done for this shot. Exiting")
            exit(0)
        else:
            with open(running_status_file, "w") as f:
                f.write("Running.\n")
                f.flush()
    except Exception as e:
        print(f"{shotid} status check, Exception: {e}\n\n{traceback.format_exc()}")

    # this is the hetdex wavelength grid
    def_wave = np.arange(3470, 5542, 2)

    # open the h5 file and load the full-sky models for the 3 exposures
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

    # make a dictionary with the full-sky models interpolated onto the hetdex wavelength grid
    # the key to the dictionary is the exposure number (1,2,3)
    ffskymod = {}
    for exp_num, model in enumerate([ffskymod_1, ffskymod_2, ffskymod_3]):
        wave, spec = np.transpose(model) # spec is in counts
        ffskymod[exp_num + 1] = np.interp(def_wave, wave, spec)
    del ffskymod_1, ffskymod_2, ffskymod_3

    fibers = get_fibers_table(shotid,survey=HDRVersion)

    # get all unique multiframes and exposures
    multiframes = np.unique(fibers['multiframe'])
    exps = np.unique(fibers['expnum'])

    if True:
        # get bad amps
        # open bad amp file and save as a dictionary like {multiframe : flag}
        config = HDRconfig(survey=HDRVersion)
        badamp_tab = fits.open(config.badamp)
        badamp_tab = Table(badamp_tab[1].data)
        # columns: 'shotid', 'multiframe', 'flag'
        badamp_tab = badamp_tab[badamp_tab['shotid'] == shotid]
        badamp_dict = {badamp_tab['multiframe'][i]: badamp_tab['flag'][i] for i in range(len(badamp_tab))}


        # save relative residual for each multiframe in each exposure
        #other than a problem with badamp_dict, any exception here is fatal, just let it abend
        for expnum in exps:
            for multiframe in multiframes:

                # ignore if flagged as bad.
                try:
                    if (multiframe in badamp_dict.keys()) and (badamp_dict[multiframe] == 0):
                        #record as bad (skipped on purpose)
                        save_filename = save_filename_pattern.format(shotid, multiframe, expnum)
                        with open(f"{save_filename}.bad","w") as ff:
                            ff.write("bad amp\n")
                        continue
                except:
                    #likely not in the badamp_dict
                    pass

                # 'here' marks all fibers belonging to this multiframe and exposure
                here = (fibers['expnum']==expnum) * (fibers['multiframe']==multiframe)

                # this is the factor to convert the calibrated spectrum to counts
                # because the full-sky model is given in counts
                gain = fibers['calfib'][here] / fibers['calfib_counts'][here]

                # calculate the relative residual
                # ratio of calfib_ffsky in counts divided by full-sky model
                res_here = fibers['calfib_ffsky'][here] / gain / ffskymod[expnum][None,:]

                # save as h5 file.
                save_filename = save_filename_pattern.format(shotid, multiframe, expnum)
                with h5py.File(save_filename, "w") as ff:
                    ff['rescor'] = res_here


        #assume all is well
        try:
            if os.path.isfile(running_status_file):
                os.remove(running_status_file)
        except Exception as e:
            print(f"{shotid} Remove .running file, Exception: {e}\n\n{traceback.format_exc()}")

        try:
            with open(done_status_file,"w") as f:
                f.write("Done.\n")
                f.flush()
        except Exception as e:
            print(f"{shotid} Create .done file, Exception: {e}\n\n{traceback.format_exc()}")
