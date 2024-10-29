# run example:
# python get_improved_ffskysub.py 20220103007 1 0

################################
#from Maja Niemeyer
# minimal changes DDavis to work in new environment
################################

import numpy as np
import sys
import os
import h5py
import traceback
import tables as tb
from astropy.table import Table

from hetdex_api.shot import open_shot_file, get_fibers_table

HDRVersion = "hdr5"

basepath="./output/ffskysub/"



def read_rescor(month, multiframe, nomask=False):

    """
    This reads the rescor file if it exists
    or replaces it with zeros if it doesn't.
    """

    filename = os.path.join(basepath,"rescor")
    if nomask:
        filename = os.path.join(filename, 'nomask')
    else:
        filename = os.path.join(filename, "mask")
    filename = os.path.join(filename, f"{month}/combined/combined_{multiframe}.h5")

    if os.path.exists(filename):
        with h5py.File(filename) as ff:
            rescor = ff['rescor'][:]
    else:
        print(f"Could not find file {filename}. Replacing with zeros.")
        rescor = np.zeros((112, 1036))
    return rescor



if __name__ == '__main__':

    # input #1 is the shotid
    # input #2 is 'nomask', i.e. only do badamp masking (recommended)
    # input #3 is whether to save it with h5py (1) or pytables (0, recommended)
    shotid = int(sys.argv[1])
    month = str(shotid)[:-5]
    nomask = int(sys.argv[2]) == 1
    save_h5py = int(sys.argv[3]) == 1

    # get output directory and filename.
#    save_filename = "maja_n/ffskysub/improved_spectra/"
    save_filename = os.path.join(basepath,"improved_spectra")

    if not os.path.isdir(save_filename):
        os.mkdir(save_filename)


    if not save_h5py:
        save_filename = os.path.join(save_filename, "pytables")
        if not os.path.isdir(save_filename):
            os.mkdir(save_filename)

    if nomask:
        save_filename = os.path.join(save_filename, "nomask")
    else:
        save_filename = os.path.join(save_filename, "mask")

    if not os.path.isdir(save_filename):
        os.mkdir(save_filename)


    #update actual file name: rcdatevshot.h5
    fn = f"rc{str(shotid)[0:8]}v{str(shotid)[8:]}.h5"
    save_filename = os.path.join(save_filename, fn)
    print(f"save_filename={save_filename}")
    
    # define hetdex wavelength grid
    def_wave = np.arange(3470, 5542, 2)

    # read full-sky model from h5 file for each exposure
    fileh = open_shot_file(shotid,survey=HDRVersion)
    try:
        ffskymod_1 = fileh.root.FullSkyModel.exp01.read()
        ffskymod_2 = fileh.root.FullSkyModel.exp02.read()
        ffskymod_3 = fileh.root.FullSkyModel.exp03.read()
    except Exception as E:
        print(f"Exception opening exposures for {shotid}. Fatal.: {E}\n\n{traceback.format_exc()}")
        try:
            if os.path.exists(fileh.filename):
                print(f"Size {shotid}, {fileh.filename}: {os.stat(fileh.filename).st_size / (1024 ** 3):0.2f} GB")
            else:
                print(f"Could not find file {fileh.filename}.")
        except:
            pass
        exit(-1)
    fileh.close()
    # interpolate full-sky model onto hetdex wavelength grid and define dictionary
    # keys are exposure numbers (1,2,3)
    ffskymod = {}
    for exp_num, model in enumerate([ffskymod_1, ffskymod_2, ffskymod_3]):
        wave, spec = np.transpose(model) # spec is in counts
        ffskymod[exp_num + 1] = np.interp(def_wave, wave, spec)
    del ffskymod_1, ffskymod_2, ffskymod_3

    # open the fiber data
    try:
        fibers = get_fibers_table(shotid, survey=HDRVersion)
    except Exception as E:
        print(f"Exception calling get_fibers_table({shotid}). Fatal.: {E}\n\n{traceback.format_exc()}")
        exit(-1)

    multiframes = np.unique(fibers['multiframe'])
    exps = np.unique(fibers['expnum'])
    fiber_ids = fibers['fiber_id']


    if len(exps) != 3:
        print(f"Error! Unexpected number of exposures in {shotid}. {len(exps)} exposures as {exps}. Fatal.")
        exit(-1)

    # calculate conversion from calibrated spectra to counts
    gain = fibers['calfib'] / fibers['calfib_counts']

    # for each multiframe, calculate the residual correction
    # which is rescor * full-sky model in calibrated units
    rescors = np.zeros(fibers['calfib_ffsky'].shape)
    for multiframe in multiframes:
        rescor = read_rescor(month, multiframe, nomask=nomask)
        for expnum in exps:
            here = (fibers['expnum']==expnum) * (fibers['multiframe']==multiframe)
            if np.sum(here) == 0:
                continue
            try:
                rescors[here] = rescor * gain[here] * ffskymod[expnum][None,:]
                print(f"Done with {multiframe} exp{expnum}.")
            except Exception as E:
                print(f"Exception computing rescor for {shotid}_{multiframe}_{expnum} : {E}\n\n{traceback.format_exc()}")

    print("Remaining zeros: {}".format(np.sum(rescors==0)))

    # subtract the residual correction from the calfib_ffsky
    improved_ffsky = fibers['calfib_ffsky'] - rescors

    # if save_h5py, save as array in h5 file
    if save_h5py:
        with h5py.File(save_filename, "w") as ff:
            ff['calfib_ffsky_rescor'] = improved_ffsky
        
    # if not, save with pytables
    # as a table with fiber_id and calfib_ffsky_rescor
    else:
        fileh = tb.open_file(
            save_filename,
            "w",
            title="residual-corrected ffsky fibers",
        )
        flags = fileh.create_table(
            fileh.root,
            "calfib_ffsky_rescor",
            Table(
                [fiber_ids, improved_ffsky],
                names=["fiber_id", "calfib_ffsky_rescor"],
            ).as_array(),
        )
        flags.flush()

        flags.cols.fiber_id.create_csindex()
        flags.flush()
        fileh.close()

    print(f"Saved to {save_filename}.")
