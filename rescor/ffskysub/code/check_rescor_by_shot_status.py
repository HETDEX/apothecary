"""

scan the multiframes (IFU+AMP) to see that they exist for all multiframes and exposures and report status if
any are missing

"""
#from datetime import datetime #, timedelta
#time_start = datetime.now()

import numpy as np
import sys
import os

import tables
from hetdex_api.survey import FiberIndex
from hetdex_api.config import HDRconfig

from astropy.table import Table
from astropy.io import fits

import traceback

from tqdm import tqdm



hdr_name = "hdr4"

check_badamps = True #do not mark as missing if the amp is in the badamps dictionary

# the input is the shotid

args = list(map(str.lower,sys.argv)) #python3 map is no longer a list, so need to cast here

#there is a list of shots to check
all_shots = []
if ("-f" in args):
    i = args.index("-f")
    if i != -1:
        try:
            all_shots = np.loadtxt(sys.argv[i + 1],dtype=int)
        except:
            print("Exception processing command line for -f")
            exit(-1)

else:
    all_shots = [int(sys.argv[1])]


shoth5_dir = f"/scratch/projects/hetdex/{hdr_name}/reduction/data"
#working_dir = "/scratch/03261/polonius/rescor/code/rescor/ffskysub/code"
working_dir = "./"
# make sure that the directory to save files to exists
#    save_filename_dir = f"maja_n/ffskysub/rescor/nomask/{month}/individual/"


#load the fiber index for this shot ... it is fast and has the minimum info we need
try: #this may fail for newer shots that are not yet in an index
    FibIndex = FiberIndex()

    for shotid in tqdm(all_shots):

        q_shot = shotid
        month = str(shotid)[:-5]

        save_filename_dir = os.path.join(working_dir, f"output/ffskysub/rescor/nomask/{month}/individual/")

        missing_status_file = os.path.join(save_filename_dir, f"rescor_by_shot_nomask_{shotid}.missing")
        done_status_file = os.path.join(save_filename_dir, f"rescor_by_shot_nomask_{shotid}.done")

        multiframes = np.unique(FibIndex.hdfile.root.FiberIndex.read_where("shotid==q_shot",field="multiframe")).astype(str)
        exps = [1, 2, 3]
        missing_mf = []

        if len(multiframes)==0: #there is a problem ... this shot was not found
            # try to load from a shot h5 file?
            try:
                dvs = str(shotid)[0:8]+"v"+str(shotid)[-3:]+".h5"
                h5 = tables.open_file(os.path.join(shoth5_dir,dvs))
                multiframes = np.unique(h5.root.Data.FiberIndex.read(field="multiframe")).astype(str)
                h5.close()
            except Exception as e:
                print(f"{shotid} status check, Exception: {e}\n\n{traceback.format_exc()}")
                # assume missing
                missing_mf.append(f"ALL_{shotid}")

        if len(multiframes) > 0 and check_badamps:
            config = HDRconfig()
            badamp_tab = fits.open(config.badamp)
            badamp_tab = Table(badamp_tab[1].data)
            # columns: 'shotid', 'multiframe', 'flag'
            badamp_tab = badamp_tab[badamp_tab['shotid'] == shotid]
            badamp_dict = {badamp_tab['multiframe'][i]: badamp_tab['flag'][i] for i in range(len(badamp_tab))}
        else:
            badamp_dict = {}

        for mf in multiframes: #walk down the list of multiframes and make sure they all exit

            if (mf in badamp_dict.keys()) and (badamp_dict[mf] == 0):
                continue #skipping since is in badamps list

            for exp in exps:
                check_file = os.path.join(save_filename_dir,f"{shotid}_{mf}_{exp}.h5")
                try:
                    if os.path.isfile(check_file) and os.stat(check_file).st_size > 10:
                        continue
                    else: #this one is missing
                        missing_mf.append(os.path.basename(check_file))
                except Exception as e:
                    print(f"{shotid} status check, Exception: {e}\n\n{traceback.format_exc()}")
                    #assume missing
                    missing_mf.append(os.path.basename(check_file))

        if len(missing_mf) > 0:
            print(f"{shotid} missing {len(missing_mf)}.")
            with open(missing_status_file,"w") as f:
                for mf in missing_mf:
                    f.write(f"{mf}\n")
                f.flush()
        else:
            #print(f"{shotid} all done.")
            with open(done_status_file, "w") as f:
                f.write("Done.\n")
                f.flush()

        #time_stop = datetime.now()

        #print(f"{shotid} time: {time_stop-time_start}")
except Exception as e:
    print(f"{shotid} status check, Exception: {e}\n\n{traceback.format_exc()}")