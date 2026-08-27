"""
Perform science reduction (astrometry, flux calibration, line and continuum detection, elixer )
  on a single shot (field). Can be HETDEX dither-style or not.

  Copies code source and works in the current directory

  input: [required] shotid (or datevshot)
         [optional] -clean  (clean up workfiles, leaving only the output and logs)
         [optional] -overwrite (delete and overwrite the datevshot directory)
         [optional] -exp <##> (specify a single exposure to reduce, if there is more than one in the shot; 0 = all)

!!! Warning !!! Not recommended for use with Pipes (have seen issues with lower level scripting calls). That is,
  when running from an idev session DO NOT to something like: python reduce_shot.py (XXXX params) | tee out.txt

This is just a large Python script. There is no defined main(), but the principle logic begins in a
  "main" commented section.

Error control (at least for now) is deliberately limited as I want no hidden errors. Mostly if anything is wrong
   I want this to break.

"""

# start limited tracking version ... starting from 20260513
# add check of *amp.dat avg_sky and plot on collapsed image
# tp adjustments
# 1.0.4 fix bug with swapping between GAIA, SDSS, PanStarrs on failure ; fix bug with tarfile path
# 1.0.5 add dust correction at 4540AA to summary.txt
# 1.0.6 add --linedet_filter options
# 1.0.7 add normed, observed EW filter for bright objects, tweak memory for rf1 (fitradecsp) based on total memory available
# 1.0.8 add color-coded pngs of up to 3x dithers and all 4 amps per IFU in analysis folder
# 1.0.9 set color-coded pngs to run with step 01 and/or as part of more complete analysis in step 04
# 1.0.10 fix issue with --hetdex tar copy sticking on the directory instead of the tar file; check lib_calib before reducing
# 1.0.11 extra status.warn conditions, IFU analysis adjustments based on exptime
# 1.0.12 add options for --linedet_parms and --force
# 1.0.13 updates to reload warn and fail conditions after a late --resume
# 1.0.14 restrict minimum SNR of *.mc raw line detections that make it into HETDEX_API based on avg_sky
# 1.0.15 fix issue where progress.dat could be overwritten as reset back to start; all progress lost
# 1.0.16 add --repair switch

__version__ = '1.0.16'

import numpy as np
import sys
import os
import glob
import shutil
from pathlib import Path
from dataclasses import dataclass
import tarfile as tar
import json
from datetime import datetime, timedelta
import multiprocessing
import subprocess
from tqdm import tqdm

try:
    from filelock import FileLock
except:
    print("You need to install filelock (e.g.: pip install --user filelock) ")
    exit(-1)

import tables
from astropy.table import Table, join, Column, MaskedColumn # unique, vstack, hstack
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord

# noinspection PyUnresolvedReferences
from h5tools import amp_stats as AmpStats
# noinspection PyUnresolvedReferences
import hetdex_tools.fof_kdtree as fof
from hetdex_api.extinction import deredden_spectra

# noinspection PyUnresolvedReferences
from elixer import global_config as G
# noinspection PyUnresolvedReferences
from elixer import spectrum_utilities as SU
# noinspection PyUnresolvedReferences
from elixer import utilities as Utils

from PIL import Image, ImageDraw, ImageFont


import importlib.util

import traceback
import psutil
import time
import copy
import matplotlib
matplotlib.use('agg')

from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.style.use('default')

########################################################################
# CONFIGURATION
########################################################################

CodeArchiveDir = "/work/03261/polonius/hetdex/single_shot/"

CorralMaxCopyLimit = 10
WorkMaxCopyLimit = 20
EchoCmds = True  #if True echo system commands to the log
FilterDetsOnBadAmps = True # if True, do NOT pass detections that are on reported bad amps to elixer for processing
DefaultClean = 1 #clean 0 does nothing, 1 cleans script files and temporary stuff, 2,3,4,5 are increasingly agressive
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot
LastKnownFplane = "fp20240731"
HETDEX_API_SNR_Thresh = 4.5 #do not include line sources where the S/N < this value
GuiderFWHM_ALL = True #if True and using the GUIDER FWHM, use all within the observation timeframe
                      #if False, just use the two nearest in time to the end of the observation

#if we know the number of active shots, use these limits / # of active shots
MaxTotalProcs_mp_rcal = 80 #flux calibration
MaxTotalProcs_mp_rf1 = 65 #line detections, 4 is about right (on averaged) to avoid memory issues on lonestar6

#10 seems to be the max for a vm-small
MaxPerShotProcs_mp_rcal = 10 #even if more Total processes are available, limit to this for any given shot
MaxPerShotProcs_mp_rf1 = 10  #even if more Total processes are available, limit to this for any given shot

#with no information, use these limite
NumProcs_mp_rcal = 6 #flux calibration
NumProcs_mp_rf1 = 4 #line detections, 4 is about right (on averaged) to avoid memory issues on lonestar6


ScriptRepo = "/work/03261/polonius/hetdex/single_shot"
LocalScriptRepo = "./local_script_repo" #useful if running multiple single shots ... can copy remotely once
                                        #then copy locally from here for each shot
                                        #set to None if you do NOT want to use a local script dir cache
                                        #  and force a copy from the main repo each time

Lock_mutex_fn = "lsr.lock"  #all instaces under this directory share this
Node_basedir = "/tmp/hetssr"
Lock_tmp_mutex_fn = "/tmp/tmp_hetssr.mutex" #all instaces on one node (common /tmp) share this
WorkDirRoot = "./"
#user specific
red1path = None # if None, will use the local (cwd) as the basepath, otherwise user can edit and specify one here
                # "/scratch/03261/polonius/red1/reductions/"

#HETRaw_archive = "/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"
HETRaw_archives = ["/work/03946/hetdex/maverick","/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"]
#HETDEXSurvey = "/corral-repl/utexas/Hobby-Eberly-Telesco/hdr5/survey/survey_hdr5.h5"
HETDEXSurvey = "/scratch/projects/hetdex/hdr5/survey/survey_hdr5.h5"
HET_by_date = "/work/03946/hetdex/maverick/"
karlgettar = "/work/00115/gebhardt/maverick/gettar/"
karlfplane = "/work/00115/gebhardt/maverick/fplane/"
karlhome = "/home1/00115/gebhardt"
hetdex_projects_path = "/scratch/projects/hetdex/"
elixer_path = os.path.dirname(importlib.util.find_spec("elixer").origin)
hetdex_api_path = os.path.dirname(importlib.util.find_spec("hetdex_api").origin)
#there is an extra "hetdex_api" at the end that points into the lower level directory for that.
#h5tools is actually a sibling
hetdex_api_path = "/".join(hetdex_api_path.split("/")[0:-1])

MAX_SAFE_AVG_SKY = 500.0 #avg sky "background" from *amp.out (d*amp.dat) above this value is a problem
FAIL_AVG_SKY = 9999.0 #avg sky "background" from *amp.out (d*amp.dat) above this value is a full fail

DIAG_AMP_IMG_VMIN_VMAX = (-30, 30)  # fixed ranges for the IFU+amp diag images (shot_analysis() and make_amp_images())
DIAG_AMP_IMG_DPI = 100

DEFAULT_MIN_SNR_FOR_ELIXER = 4.3 #can be changed by --linedet_filter switch (and might not agree with HETDEX_API_SNR_Thresh)
#pretty generous here, maybe should also adjust the warnings based on exptime and observing conditions?
WARN_NUM_LINE_DETS = 5000
WARN_NUM_CONT_DETS = 500

#since this is done outside (that is, the SLURM sets this up) it is not useful here, in this code, to know this
# A way to deal with it here, would be to change the mutex to a sempahore and have a resource count
#   and then throttle any use, regardless of SLURM ... e.g for a given /tmp/ only allow 1 shot to run for vm-small
#   or 10-12 or so for normal .. any additionals would wait on the semaphore
AssumedMemFootprint = 20.0  # GB assumption
SafeActiveShotsSleep = 30.0 # recheck every xx seconds
try:
    ApproxBaseRAM = psutil.virtual_memory()[0] / (1024**3) #in GB (e.g. ~32GB for vm small, 256 GB for normal on LS6)
    #need somewhere around 20GB for normal big shots (once IFU is full)
    #varies depending also on the number of exposures, but 20GB is a safe rule of thumb
    MaxSafeActiveShots = int(ApproxBaseRAM//AssumedMemFootprint)
    print(f"*** setting MaxSafeActiveShots to {MaxSafeActiveShots}. "
          f"BaseRAM {ApproxBaseRAM:0.1f}GB, Footprint ~ {AssumedMemFootprint}GB")


    #make RAM based adjustments ... assumption is 251.1GB available for development or normal
    #but if vm-small, this is only 31.2
    mem_mux = max(1.0,round(251.1 / ApproxBaseRAM)) #round to an integer, but cast as float and at least 1

    MaxTotalProcs_mp_rcal = max(1,round(MaxTotalProcs_mp_rcal / mem_mux))  # flux calibration
    MaxTotalProcs_mp_rf1 = max(1,round(MaxTotalProcs_mp_rf1 / mem_mux)) # line detections

    # MaxPerShotProcs_mp_rcal = max(1,round(MaxPerShotProcs_mp_rcal / mem_mux))
    # MaxPerShotProcs_mp_rf1 = max(1,round(MaxPerShotProcs_mp_rf1 /mem_mux))

except:
    ApproxBaseRAM = -1
    MaxSafeActiveShots = 0


#execute steps
s01_run1s = True
s01b_amp_images = True

s02_vdrp = True
do_panstarrs = False #only run PanSTARRS if true, otherwise just run the usual GAIA and SDSS

s03_fluxcal = True

s04_make_shot = True
s04a_get_ifucens = True
s04b_rfft = True  #rfft = run full field on tmp
s04c_rcal_all = True
s04d_shot_h5 = True
s04e_amp_stats = True  #includes updating FiberIndex, FiberMask and AmpStats
s04f_analysis = True
s04_make_shot = s04_make_shot | s04a_get_ifucens | s04b_rfft | s04c_rcal_all | s04d_shot_h5 | s04e_amp_stats | s04f_analysis#sanity catch

s05_detection = True
s05b_rdet_rf1 = True
s05c_rgetmax = True
# s05d_detection_tables = False  #builds as fits files... this is replaced by the hdf5 version
s05e_detection_hdf5 = True
s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05e_detection_hdf5 #sanity catch
#s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05d_detection_tables | s05e_detection_hdf5 #sanity catch

s06_catalogs = True
s06b_fof = True      #cluster the lines and continuum sources (separately)
s06c_diagnose = True #run Diagnose
s06d_elixer = True   #run elixer
s06e_source_cat = True #make a source catalog

s06_catalogs = s06_catalogs | s06b_fof | s06c_diagnose | s06d_elixer | s06e_source_cat

#
# if False: #testing
#     print("#################### TESTING ##########################")
#     # execute steps
#     s01_run1s = False
#
#     s02_vdrp = False
#     do_panstarrs = False  # only run PanSTARRS if true, otherwise just run the usual GAIA and SDSS
#
#     s03_fluxcal = False
#
#     s04_make_shot = False
#     s04a_get_ifucens = False
#     s04b_rfft = False
#     s04c_rcal_all = False
#     s04d_shot_h5 = False
#     s04e_amp_stats = False
#     s04_make_shot = s04_make_shot | s04a_get_ifucens | s04b_rfft | s04c_rcal_all | s04d_shot_h5 | s04e_amp_stats  # sanity catch
#
#     s05_detection = True
#     s05b_rdet_rf1 = False
#     s05c_rgetmax = True
#     #s05d_detection_tables = False  #builds as fits files... this is replaced by the hdf5 version
#     s05e_detection_hdf5 = True
#    # s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05d_detection_tables | s05e_detection_hdf5  # sanity catch
#     s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05e_detection_hdf5  # sanity catch
#
#     s06_catalogs = True
#     s06b_fof = True  # cluster the lines and continuum sources (separately)
#     s06c_diagnose = True  # run Diagnose
#     s06d_elixer = True  # run elixer
#     s06e_source_cat = False  # make a source catalog
#
#     s06_catalogs = s06_catalogs | s06b_fof | s06c_diagnose | s06d_elixer | s06e_source_cat





########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:
    args = None

    hetdex: bool = False #if True, can re-run hetdex shots, otherwise they are not allowed
    het_raw_path: str = "" #path to local het_raw instead of under cwd()
    clean: int = 0 #post run clean level; 0 = do not clean
    clean_only: bool = False #if True, do nothing except for -clean
    clean_done : bool = False #set to True if the post_clean() has been performed
    #simul : int = 0 #number of simultaneous shots being run (e.g. tasks per node), 0 is unset
    update_local_repo: bool = False
    update_only : bool = False
    overwrite: bool = False
    shot_only: bool = False
    resume: bool = False
    shotid: int = 0
    datevshot: str = ""
    exp: int = 0  #specific exposure number to reduce
    email: str = ""
    numexp: int = 0 #number of exposures in the shot
    total_exp_time : float = 0.0 #in seconds
    dither_exp_times = []
    cwd_orig: str = os.getcwd()
    cwd: str = os.getcwd()
    virus_tar_path : str = None
    #red1dir: str = f"{os.getcwd()}" #e.g. <path>/red1/reductions
                    # BUT "reductions" is added later, so stop at the "red1" equivalent; so is same as just cwd_orig here
    scriptdir: str = ""
    gettar_fn: str =  ""  # the runs* or runt* file from karlgettar folder with the date, shot, exp data
    starcat_ast = "gaia" #gaia, sdss, panstarrs   for Astrometry
    starcat_cal = "sdss" #gaia, sdss, panstarrs   for FluxCalibration

    orig_stdout = None
    orig_stderr = None
    file_stdout = None

    guider_fwhm = None

    shot_obj = None #a name of sorts
    special = 0 #do some special, direct edit code stuff
    hetdex_original = False #set to True if this shot is in the original hetdex data
                            ## (vs the 'hetdex' member which is true if --hetdex is specified to allow this)
    multifits_only = False  #stop once the mutli*fits files have been generated
    mf_file_status = {}
    sub_shot = None #special case useage, use this shotnum instead of that from the datevshoot for some IFU lookups
    ifuslot = None
    build_rta = False #set to true if we need an rta file before calling into vdrp
    shot_ra = -999.9
    shot_dec = -999.9
    shot_track = -1 #0 = east , 1 = west
    avg_sky = None #average sky after run_run1s values > 1000.0 are a problem
                   #note: though related too, this is NOT exactly the same as the h5 AmpStats table "avg_orig" column
    code_fn_to_copy = None
    set_warn = False
    dither_configuration = None #-1 is bad, 0 = non-standard (maybe not dithered), 1 = standard hetdex
    dither_norms = []
    dither_positions = []
    linedet_filter = 0 #default, normal filtering
    linedet_parms = [1, 0.0] #old HETDEX style = (3, 0.5) #must be (<int>,<float>)

    made_amp_images = False #used to indicate the diagnostic IFU+amp images were already made in this run instance
    #active_ifus = None #number of active IFUs (note not normally set until late into step04e?)
    amp_stats_problem = -1
    num_bad_amps = -1
    num_all_amps = -1
    ifu_list = [] #as multiframe
    ifu_linedet_ct = [] #actual count of line dets in the requested selection
    ifu_linedet_ratio = [] #ratio to what was maximally expected
    mf_clean_image_avg_row = None

    num_line_dets = -1 #unset
    num_cont_dets = -1 #unset
    ratio_line_dets = 1.0 #unset, assume normal
    snr_rescale = 1.0 #recommended multiplier to raise the normal, minimum SNR floor (based on long exposure noise)

    # for cleanup at the end
    copy_lock_file = None
    node_clean_done = False
    dtprog = None
    write_progress = True #allow progress.dat to be written out

########################################################################
# Basic user input
########################################################################


def check_version(cfg):
    """
    really basic check ... does this file's version match the "archived" version?

    :return:
    """

    try:
        fn = os.path.join(CodeArchiveDir,"reduce_shot.py")
        arc_ver = None
        if os.path.exists(fn):
            with open(fn,"r") as f:
                found_version = False
                #near the top
                while not found_version:
                    line = f.readline()
                    if not line:
                        break
                    if "__version__" in line:
                        toks = line.split('=')
                        arc_ver = str(toks[1].replace('\'', "").replace('\"', "").strip())
                        found_version = True
        else:
            print(f"Could not check version. Counld not locate code archive: {fn}")

        if arc_ver is not None:
            if arc_ver != __version__:

                cfg.code_fn_to_copy = fn
                print(f"**********************************************************************************************")
                print(f"* Warning! Local version {__version__} does not match code archive version {arc_ver} at:")
                print(f"* \t{fn}")
                print(f"* You may want to cancel and update reduce_shot.py from the above path.")
                print(f"* If you do not cancel, this code will continue as is in 30 seconds ... ")
                print(f"**********************************************************************************************")
                for _ in tqdm(range(30)):
                    time.sleep(1.0)
            # else:
            #     print(f"Confirmed. Local version {__version__} matches code archive version {arc_ver} at {fn}")

    except:
        print(f"Could not check version. Exception!",traceback.format_exc())

print(f"version {__version__}")
print(f"Command: {' '.join(sys.argv)}")

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

cfg = Config()
cfg.args = copy.copy(args) #need to check some args if present, later on

check_version(cfg)

if "-help" in args:
    #do NOTHING else except print the help
    help = """
    
    usage: python reduce_shot.py [switches] <shot>
           where <shot> is YYYYMMDDSSS or YYYYMMDDvSSS
    
    switches (all optional):
    --cal_flux <str>: Choose one of "gaia" or "sdss" as the primary star catalog for flux calibration. The other becomes
                      the secondary and PanSTARRS is always the tertiary. "sdss" is the default.
                      
    --cal_astro <str>: Choose one of "gaia" or "sdss" as the primary star catalog for astrometric calibration. 
                       The other becomes the secondary and PanSTARRS is always the tertiary. "gaia" is the default.
    
    --clean <integer> : -5 to 5; 0 performs no cleanup (e.g. for full debugging). Default = (1) if not specified.
                      negatives do the same as positives except they force the clean-up immediately and as the only action
                      0 = No clean up at all. Intermediate files and all tempory scripts reamin. Good for debugging.
                      1 = *** default *** Basic cleaning of scripts and non-informative intermediate files
                          (--resume can still work if accompanied by --update)
                      2 = removes additional intermediate files, and most debugging files, vdrp, etc
                          (--resume can still work if accompanied by --update)
                      3 = removes logging information, analysis, match_pngs, etc
                          (--resume can still work if accompanied by --update)
                      4 = removes the .mc, .spec, .list, etc files 
                          (!!! cannot use --resume at level 4 or after !!!) 
                      5 = removes EVERYTHING except the shot h5 file
            !!! Notice: elixer content is only removed with (-1) to (-5) or 5 as it is a post-run manual SLURM queue
                      
    --clear_mutex : Special operation - clears the concurrency mutex, resetting allowed concurrent active processes.
                    This takes priority over all other switches and will terminate with this action.
                                                          
    --exp <integer> : will operate on only the specified exposure (e.g. in a multi-exposure observation, can select
                      exactly one to reduce). If not present or set to (0), will use all exposures for the observation. 
                      
    --email <str> : if provided will attach to the elixer slurm job so this email address will get notifications
                    
    --force : continue the reduction regardless of select exit conditions
              + ignore excessive line detections 
              + ignore excessive bad amps
        
    --help : display this help text and exit
        
    --hetdex : if present, overrides the restriction on running existing HETDEX shots
    
    --ifuslot <str(3)> : Special operation - work with ONLY this IFUSlot
                         only used in conjuection with --multifits_only
    
    --linedet_filter : preset filter/selection for line detections to send to ELIXer. Default is (1)
                    0 = default: snr > 4.8, chi2 <= 2.5, 1.5 <= linewidth <= 16, chi2fib <= 4.5, continuum >= -3
                    1 = OFF + normal HETDEX_API.  Only snr >= 4.5 filter.
                    1.x+ = OFF + HETDEX_API snr set to this value (i.e. for values 1.1 and up, single decimal precision)
                    
    --linedet_parms : grid (int) and step (float) for line detections as <grid>,<step>
                      3x dither default is 3,0.5
                      1x dither default is 1,0.0
    
    --local_het_raw_path : if present, specifies the "local" het_raw path to use for locally copied data. This is a path
                     to which /het_raw/ will be appended.
    
    --multifits_only : if present, stops processing and terminates after the multi*fits files have been created.
                       Can be combined with --exp to limit which exposures to process.
                        
    --nolimit : if present, overrides the in-code limiting of simultaneous active shots per node
             
    --overwrite : removes the shot working directory completely and (re)starts fresh. 
              !!! Notice: --resume has priority over --overwrite
    
    --prep_compress : Special operation - creats an executable script to create the ssr_<shot>.h5 files, but does not
                                          execute that script. Takes a single integer parameter to set the maximum
                                          concurrent processes to use followed by the shots to include (wildcards
                                          allowed but the string should be in quotes).
    
    --queue_elixer : Special operation - Inserts provided shots (wildcard allowed, but should be in quotes)
                                         into SLURM queue as previously designated. Updates elixer slurm scripts,
                                         if needed, to adjust to path changes.
    
    --repair : basically --resume but with a repair/re-copy of files that would have been lost to cleaning
    
    --resume : (re)starts roughly at the last completed step (see sciXXXX/progress.dat)
           !!! Notice: This does NOT re-run steps that completed with failures, it only re-runs incomplete steps.
           !!! Notice: --resume has priority over --overwrite
           
    --shot_only : ONLY (re)build the shot h5 file. Do NOT run detections or elixer.
               
    --sub_shot <str(3)> : Special operation - substitute this shot number for that of the datevshot for the gettar lookup 
                          only used in conjunction with --multifits_only
               
    --update : removes and re-fetches the local_script_depo prior to running
               on a --resume, also updates the scripts already in the shot working directory
    
    """

    print(help)
    exit(0)

len_args = len(args)
queue_elixer = False
prep_compress = -1 #not just a boolean, use as the max simultaneous shots to process
FORCE_CONTINUE = False


if "-force" in args:
    FORCE_CONTINUE = True
    args.remove("-force")

if "-cal_flux" in args:
    i = args.index("-cal_flux")
    try:
        cfg.starcat_cal = str(args[i+1]).lower()
        if cfg.starcat_cal in ["sdss","gaia"]:
            pass #all is good
        else:
            print(f"Invalid -cal_flux specified: {args[i+1]} ; must be 'sdss' or 'gaia' ")
            exit(-1)
    except:
        print(f"Invalid -cal_flux specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-cal_flux")

if "-cal_astro" in args:
    i = args.index("-cal_astro")
    try:
        cfg.starcat_ast = str(args[i+1]).lower()
        if cfg.starcat_ast in ["sdss","gaia"]:
            pass #all is good
        else:
            print(f"Invalid -cal_astro specified: {args[i+1]} ; must be 'sdss' or 'gaia' ")
            exit(-1)
    except:
        print(f"Invalid -cal_astro specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-cal_astro")

if "-clear_mutex" in args:
    print("Hidden switch: clearing mutex, resetting allowed concurrent active processes.")
    print("This takes priority over all other switches and will terminate with this action.")
    fns = glob.glob(f"{Node_basedir}/*.sync")
    print(f"Clearing {len(fns)} active sync files.\n{[os.path.basename(x) for x in fns]}")
    os.system(f"rm {Node_basedir}/*.sync")

    mux_root = "/".join(os.getcwd().split("/")[:4]) + "/ssr_mux/corral"
    fns = glob.glob(f"{mux_root}/*.sync")
    print(f"Clearing {len(fns)} /corral sync files.\n{[os.path.basename(x) for x in fns]}")
    os.system(f"rm {mux_root}/*.sync")

    mux_root = "/".join(os.getcwd().split("/")[:4]) + "/ssr_mux/work"
    fns = glob.glob(f"{mux_root}/*.sync")
    print(f"Clearing {len(fns)} /work sync files.\n{[os.path.basename(x) for x in fns]}")
    os.system(f"rm {mux_root}/*.sync")

    print("Done. Exiting.")
    exit(0)

if "-queue_elixer" in args:
    # this is option is to be used on its own call, not part of a normal reduction
    # it looks for the shot under the current directory (reduction must already be completed in a previous run)
    #  and calls sbatch for the line and conts (again, which should have been previously created)
    # e.g. it just makes it simpler to call sbatch on already created slurm jobs
    #      it shoud be executed from the normal login node
    print("Hidden switch : queueing elixer slurm jobs that match datevshot ...")
    print("Usage: python reduce_shot.py --queue_elixer <datevshot>")
    print("Usage: wildcards allowed, but put in quotes and DO NOT prefix with 'sci ")
    args.remove("-queue_elixer")
    queue_elixer = True

if "-prep_compress" in args:
    # this is option is to be used on its own call, not part of a normal reduction
    # it looks for the shot under the current directory (reduction must already be completed in a previous run)
    #  and calls sbatch for the line and conts (again, which should have been previously created)
    # e.g. it just makes it simpler to call sbatch on already created slurm jobs
    #      it shoud be executed from the normal login node
    print("Hidden switch : preparing default SSR compression calls that match datevshot ...")
    print("Usage: python reduce_shot.py --prep_compress [max_simultaneous, default=3] <datevshot>")
    print("Usage: wildcards allowed, but put in quotes and DO NOT prefix with 'sci ")

    i = args.index("-prep_compress")
    try:
        prep_compress = int(args[i+1])
        if prep_compress < 0:
            print(f"Invalid --prep_compress specified. Must be non-negative")
            exit(-1)
    except:
        print(f"Invalid --prep_compress specified")
        exit(-1)

    del args[i + 1]
    args.remove("-prep_compress")
else:
    prep_compress = -1

if "-clean" in args:
    i = args.index("-clean")
    try:
        cfg.clean = int(args[i+1])
        if cfg.clean < 0:   #negative values are the same, but force just the clean operation (e.g. to be used on an old run)
            cfg.clean *= -1
            cfg.clean_only = True
    except:
        print(f"Invalid --clean specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-clean")
else:
    cfg.clean = DefaultClean #usually this is level 1


if "-linedet_filter" in args:
    i = args.index("-linedet_filter")
    try:
        cfg.linedet_filter = float(args[i+1])
        if cfg.linedet_filter not in [0,1]:
            if cfg.linedet_filter >= 1.1:
                pass #specific value
            else:
                print(f"Invalid --linedet_filter specified. Must be in [0,1] or >= 1.1. See --help")
                exit(-1)
    except:
        print(f"Invalid --linedet_filter specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-linedet_filter")
else:
    cfg.linedet_filter = 0 #default


if "-linedet_parms" in args:
    i = args.index("-linedet_parms")
    try:
        remove_chars = "()[]<>{} "
        toks = args[i+1].translate(str.maketrans("", "", remove_chars)).split(",")
        if len(toks) != 2:
            print(f"Invalid --linedet_filter specified: {args[i+1]}")
            exit(-1)

        cfg.linedet_parms = (int(toks[0]),float(toks[1]))

    except:
        print(f"Invalid --linedet_parms specified {args[i+1]}")#,traceback.format_exc())
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-linedet_parms")


# if "-simul" in args:
#     i = args.index("-simul")
#     try:
#         cfg.simul = int(args[i + 1])
#         if cfg.simul < 0:
#             print(f"Invalid -simul specified")
#             exit(-1)
#     except:
#         print(f"Invalid -simul specified")
#         exit(-1)
#
#     del args[i + 1]  # args.pop(0) #remove THIS file
#     args.remove("-simul")
# else:
#     cfg.simul = 0  # usually this is level 1

if "-update" in args:
    cfg.update_local_repo = True
    args.remove("-update")

if "-overwrite" in args:
    cfg.overwrite = True
    args.remove("-overwrite")

if "-shot_only" in args:
    cfg.shot_only = True
    args.remove("-shot_only")
    #note: reminder, if you want the multi*fits, they are in the ../reductions directory for temporary holding

if "-multifits_only" in args:
    cfg.multifits_only = True
    args.remove("-multifits_only")



if "-exp" in args:
    i = args.index("-exp")
    try:
        cfg.exp = int(args[i+1])
    except:
        print(f"Invalid -exp specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-exp")


if "-repair" in args: #like resume, but with a copy of files
    cfg.repair = True
    args.remove("-repair")

if "-resume" in args: #opposide of --overwite ... do NOT touch the (intermediate) output of the working directory
    cfg.resume = True
    args.remove("-resume")

if "-exp" in args:
    i = args.index("-exp")
    try:
        cfg.exp = int(args[i+1])
    except:
        print(f"Invalid -exp specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-exp")

if "-email" in args:
    i = args.index("-email")
    try:
        cfg.email = args[i+1]
        #really basic sanity check
        if '@' not in cfg.email:
            print(f"Invalid -email format specified: {cfg.email}")
            exit(-1)
    except:
        print(f"Invalid -email specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-email")


if "-hetdex" in args:
    cfg.hetdex = True
    args.remove("-hetdex")

if "-local_het_raw_path" in args:
    i = args.index("-local_het_raw_path")
    try:
        cfg.local_het_raw_path = args[i+1]
        del args[i + 1]
        args.remove("-local_het_raw_path")

        if not os.path.exists(cfg.local_het_raw_path):
            Path(cfg.local_het_raw_path).mkdir(parents=True, exist_ok=True)

    except:
        print(f"Invalid local_het_raw_path path specified: {args[i+1]}")
        exit(-1)
else:
    cfg.local_het_raw_path = None


if "-nolimit" in args:
    MaxSafeActiveShots = 0
    print("*** --nolimit Override. Do NOT enforce simultaneous active shot limit per node.")
    args.remove("-nolimit")


if "-special" in args:
    i = args.index("-special")
    cfg.special = int(args[i+1])
    del args[i + 1]
    args.remove("-special")
    print(f"*** --special condition invoked. Condition = {cfg.special}")

if "-sub_shot" in args:
    i = args.index("-sub_shot")
    try:
        cfg.sub_shot = args[i+1]
        if len(cfg.sub_shot) !=3 and 1 <= int(cfg.sub_shot) <= 999:
            print(f"Invalid -sub_shot specified: {args[i + 1]}")
            exit(-1)
    except:
        print(f"Invalid -sub_shot specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-sub_shot")

if "-ifuslot" in args:
    i = args.index("-ifuslot")
    try:
        cfg.ifuslot = args[i+1]
        if len(cfg.ifuslot) !=3 and 1 <= int(cfg.ifuslot) <= 999:
            print(f"Invalid -ifuslot specified: {args[i + 1]}")
            exit(-1)
    except:
        print(f"Invalid -ifuslot specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-ifuslot")

#whatever is left should be the shot
if len(args) != 1:
    #could be just an update or a help (but help is handled earlier)
    if cfg.update_local_repo and len_args == 1: #was just an update
        cfg.update_only = True
    elif prep_compress >= 0: #use all the sci* directores under the cwd
        cfg.shotid = None
        cfg.datevshot = None
    else:
        print(f"Fatal: Problem with remaining args: {args}")
        if "-clear" in args:
            print(f"--clear found ... did you mean --clean ?")
        print(f"exiting....")
        exit(-1)
else:
    if not (queue_elixer or prep_compress > -1):
        try:
            #might have 'v' or 'd' or 's' as the separator between date and shot number
            cfg.shotid = int(args[0].replace("v","").replace("s","").replace("d",""))
        except:
            pass

        if not (20170101000 < cfg.shotid < FutureShotDateLimit):
            print(f"Fatal: Invalid shotid: {args[0]}")
            exit(-1)

        cfg.datevshot = str(cfg.shotid)[0:8] + "v" + str(cfg.shotid)[8:]
    else:
        cfg.datevshot = args[0] #as is ... wildcards and all

#check datevshot ... can only be numeric or in datevshot format, but could be truncated

print(f"[{cfg.datevshot}] Evaluating with datevshot = {cfg.datevshot}",flush=True)

########################################################################
# worker functions
########################################################################

def check_lib_calib(cfg):
    """
    if the lib_calib/yyyymm directory is in the process of updating (e.g. if calibrations are being made)
      do not allow a reduction to start

    NOTE: this is incomplete protection as a calibration could start AFTER a reduction and it would not be
    protected.

    :param yyyymm:
    :return:
    """

    #todo: check the path to lib_calib/yyyymm and see if the update in progress marker is present
    #todo: is this even worth doing?? ... might be more meaningful since the time window is actually pretty big once
    #      you include the --multifts_only run for the whole month to build up the new reschi stuff
    rc = 0
    try:
        month = cfg.datevshot[0:6]
        path_check = os.path.join(hetdex_projects_path,f"lib_calib/{month}/status.updating")
        if os.path.exists(path_check):
            rc = 1 #warn, we are updating
    except:
        print(f"[{cfg.datevshot}] Could not check on lib_calib status.")
        rc = -1

    return rc



def safe_cd(path):
    """
    attempt to cd and return True only if new cwd() matches

    this fails for ".." for example
    so is only useful, in its current form for non wildcard or special paths

    :param path:
    :return:
    """

    try:
        orig = os.getcwd()
        os.chdir(path)
        if os.getcwd()[-1*len(path):] == path:
            return True
        else:
            os.chdir(orig)
            print(f"!!!!! WARNING. Failed to cd to {path} ")
            return False
    except:
        return False


def run_queue_elixer(cfg):
    """

    :param cfg:
    :return:
    """

    cwd = os.getcwd()
    fns = glob.glob(f"sci{cfg.datevshot}")
    print(f"Attempting to SLURM queue using pattern: {cfg.datevshot}")
    print(f"  matches: {fns}")
    for fn in fns:
        try:
            for dettype in ["out","line","cont"]:
                datevshot = os.path.basename(fn)[3:] #e.g. sci20260421v011, strip off the sci
                slurm_path = os.path.join(fn, f"elixer/{dettype}")
                if safe_cd(slurm_path):

                    # do not run if already "done"
                    if os.path.exists("elixer.done") or os.path.exists(f"elixer_{datevshot}_cat.h5"):
                        print(f"[{datevshot}] ELiXer already done.")

                    else: #if os.path.exists("elixer.slurm"):
                        if  os.path.exists("elixer.slurm") and os.path.exists("elixer.run"):
                            #check the paths and fix as needed
                            #e.g. from
                            # --shot_h5 /scratch/03261/polonius/parallel/shela/sci20240929v020/20240929v020.h5
                            # to
                            # --shot_h5 /scratch/03261/polonius/parallel/shela/set_xab/sci20240929v020/20240929v020.h5
                            # so, check the paths (--shot_h5 , --diagnose) specifically
                            # -OR-
                            # get the base path for --shot_h5 compare that to the current workind dir and
                            #  if they do not match, string repace all instanaces of the existing old basepath
                            #    with the current working dir
                            #      as there will be multiple lines and multiple paths in each line of the .run file
                            #
                            #  NOT exactly, be smarter ... probaby over kill, but may need to do this line by line
                            #     in case this is a combination of run files? And/or there are different include paths
                            #       different lines
                            #  1. check --shot_h5, --diagnose (and any other)
                            #      IF the file does not exist, reaplce the base path with the current basebase
                            #        up to "elixer/out" and try again
                            #  2. if the NEW path DOES exist, we are good, replace and move on
                            with open("elixer.run","r") as f:
                                #todo: open another file to receive the write? then rename??
                                # may be cleaner than trying to overwrite in place, since this is line by line
                                # NOTE: if not line by line, could run a global replace with AWK as I do elsewhere

                                # todo: test ... make a copy of elixer.run as a backup
                                #  modify this code to NOT actually issue the slurm command below, just print
                                #     that it would be exectued
                                #  then run this as a test and compare the new elixer.run to the saved one
                                #
                                # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits")
                                # cd dispatch_0001 ; /usr/bin/python /work/03261/polonius/maverick/science/sciscripts/elixer.test/elixer/elixer.py  -f --slurm 0 --nodes 1 --log info --shot_h5 /scratch/03261/polonius/parallel/shela/sci20241005v013/20241005v013.h5 --diagnose /scratch/03261/polonius/parallel/shela/sci20241005v013/diagnose_classifications.tab --png --error 3.0 --neighborhood 10.0 --ntasks_per_node 25 --timex 0.8 --post_merge 2 --merge_name elixer_20241005v013_cat.h5 --email dustin@astro.as.utexas.edu --name out --dets line.dets --hdf5 /scratch/03261/polonius/parallel/shela/sci20241005v013/20241005v013_line.h5 --dispatch dispatch_0001 -f ; cd ..
                                #

                                #get the first row and assume everything else is much the same
                                line = f.readline()
                                line = line.split(";")[1]
                                toks = line.split()
                                idx = toks.index("--shot_h5")
                                idx += 1
                                if not os.path.exists(toks[idx]):
                                    #assume we have changed directories
                                    old_path = os.path.dirname(toks[idx])
                                    shot_fn = os.path.basename(toks[idx])
                                    leaf_dir = os.path.basename(os.path.dirname(toks[idx]))

                                    new_path = os.path.join(cwd,leaf_dir,shot_fn)
                                    if os.path.exists(new_path):
                                        #do the update
                                        new_path = os.path.dirname(new_path) #strip off the file at the end
                                        print(f"[{shot_fn.replace('.h5','')}] changing to new path ....")
                                        system_command(cfg,"cp elixer.run elixer.run.original_path")
                                        system_command(cfg, f"sed -i s#{old_path}#{new_path}#g elixer.run")

                        #don't need the cd commands anymore since using safe_cd() above
                        #cmd = f"cd {slurm_path} ; sbatch elixer.slurm ; cd {cwd}"
                        cmd = f"sbatch elixer.slurm"
                        #print(f"Sending cmd to shell: {cmd}")
                        #print(f"TESTING DUMMY: >> {cmd}")
                        system_command(cfg,cmd)
                    #end else

                    #return to previous directory
                    os.chdir(cwd)
        except:
            print(traceback.format_exc())

    os.chdir(cwd)

def run_prep_compress(cfg,max_simultaneous=3):
    """

    :param cfg:
    :return:
    """
    bg = max_simultaneous #run in the background
    cwd = os.getcwd()
    if cfg.datevshot is not None:
        print(f"Attempting to prepare default SSR compression using pattern: sci{cfg.datevshot}")
        fns = glob.glob(f"sci{cfg.datevshot}")
    else:
        print(f"Attempting to prepare default SSR compression using pattern: sci*")
        fns = glob.glob(f"sci*")

    #must be of form sciYYYYMMDDvSSS ... exclude other matches
    print(f"  matches: {fns}")

    for cpath in fns: #fn is the full path to the sciXXX directory, not the files
        try:
            datevshot = os.path.basename(cpath).replace("sci","")

            if shutil.which("ssr_hdf5"):
                cmd = f"ssr_hdf5 --compression 2"
            else:
                cmd = f"python ssr_hdf5.py --compression 2"

            if bg >= 0:
                cmd += f" --bg {bg}"

            fn = os.path.join(cpath,f"{datevshot}.h5")
            if os.path.exists(fn):
                cmd += f" --shot_h5 {fn}"
            else: # we are done, cannot build out this one
                print(f"[{fn}] not found.")
                continue

            #elixer h5
            elixer = True
            fn = os.path.join(cpath, f"elixer/out/elixer_{datevshot}_cat.h5")
            if os.path.exists(fn):
                cmd += f" --elixer_h5 {fn}"
            else: # we can run without elixer, but should warn
                elixer = False
                print(f"[{datevshot}] !WARNING! No correspondig elixer_*_cat.h5 file found. Not fatal, but "
                      f"there will be no ELiXer data or reports in the output.")

            if elixer:
                fn = os.path.join(cpath, f"elixer/out/all_pngs")
                if os.path.exists(fn):
                    cmd += f" --images {fn}"
                else:  # we can run without elixer report images, but should warn
                    print(f"[{datevshot}] !WARNING! No correspondig report images found. Not fatal, but "
                          f"there will be no ELiXer reports in the output.")

            if bg:
                cmd += " &"
            with open("compress.run","a") as f:
                print("compress.run << ",cmd,flush=True)
                cmd += "\n"
                f.write(cmd)
        except:
            print(traceback.format_exc())

    os.chdir(cwd)

def hetdex_dither(cfg):
    """
    check the dither positions if 3 exposures ... is it HETDEX ?
    since this uses the shot h5 file we generate, need to wait until that is ready to run
    (can use this to see if we should bother with detections)
    (e.g. if just one exposure OR is 3 exposures as HETDEX dither, then, yes, run detections)
    :param cfg:
    :return: -1 (error) 0 = no  1 = yes
    """

    rc  = 0
    try:
        h5 = tables.open_file(os.path.join(cfg.cwd,f"{cfg.datevshot}.h5"),mode='r')
        x = h5.root.Shot.read(field="xditherpos")[0]
        y = h5.root.Shot.read(field="yditherpos")[0]

        try:
            logstr = f"[{cfg.datevshot}] dithers: ["
            cfg.dither_positions = []
            for xx,yy in zip(x,y):
                logstr += f" ({xx},{yy})"
                cfg.dither_positions.append((xx,yy))
            logstr += " ]"
            print(logstr)
        except:
            pass

        #very roughly ... allowing for some pretty big errors to even try
        if len(x) == 3 and len(y) == 3:
            if x[0] == 0 and y[0] == 0:
                # if 0.0 < x[1] < 1.7 and -1.2 < y[1] < 0.0:
                #     if 0.0 < x[2] < 1.7 and 0.0 < y[2] < 1.5:
                #         rc = 1


                #expect 1 move in x of about 1.5" and then x2 close to x1
                #and 2 moves in y with opposite signs and separation around 1.5"
                # keep this pretty loose

                if abs(x[2]-x[1]) > abs(y[2]-y[1]):
                    #flip x and y, just for this logic
                    #it does not really matter which holds close to constant and which moves around 1.5"
                    #for the fill
                    y = h5.root.Shot.read(field="xditherpos")[0]
                    x = h5.root.Shot.read(field="yditherpos")[0]

                if 0.5 < abs(x[1]) < 2.5 and abs(x[2] - x[1]) < 0.5:
                    if y[1] * y[2] < 0 and 0.5 < abs(y[2] - y[1]) < 2.5:
                        rc = 1

                # #keep this pretty loose
                # if (0.0 < x[1] < 2.0) and (-2.0 < y[1] < 0.0) and \
                #    (0.0 < x[2] < 2.0) and ( 0.0 < y[2] < 2.0):
                #     rc = 1
    except:
        rc = -1
        print(f"[{cfg.datevshot}] Exception in hetdex_dither()", traceback.format_exc())

    try:
        h5.close()
    except:
        pass

    return rc

def post_clean(cfg):
    """
    clean up after the run ...

    !!! this is not the safest way to do this, but it is fast and simple (using system commands to rm and unlink)

    for now, there are only two differences ...
    1 = basic clean (remove script files, intermediate files)
        (keep: mc files, spec files, etc)
        (discard vdrp ?)
    2 = remove same as (1) but also more agressively remove (logs, helper data files, etc)

    :param cfg:
    :return:
    """

    try:

        if cfg.clean_done:
            return

        #clean up the copy lock file
        if cfg.copy_lock_file is not None:
            try:
                lock = FileLock(cfg.copy_lock_file)
                with lock:
                    os.remove(cfg.copy_lock_file)
                    cfg.copy_lock_file = None
            except:
                pass

        #always try to clean up /tmp
        node_clean(cfg)

        if cfg.clean <=0:
            print(f"[{cfg.datevshot}] No -clean")
            cfg.clean_done = True
            return
        else:
            print(f"[{cfg.datevshot}] --clean {cfg.clean}")


        if cfg.clean_only:  #this is at the top and cfg.cwd has NOT been set to cfg.datevshot
            if not safe_cd(os.path.join(cfg.cwd,f"sci{cfg.datevshot}")):
                print(f"[{cfg.datevshot}] Could not initiate --clean")
                return
            else:
                #we did successfully jump to sciXXXX direcotry so update to that for the remainder
                cfg.cwd = os.getcwd()
        else:
            if not safe_cd(cfg.cwd):
                print(f"[{cfg.datevshot}] Could not initiate --clean")
                return

        #os.chdir(os.path.join(cfg.cwd,f"sci{cfg.datevshot}"))
        cfg.cwd = os.getcwd()  # now under the sci<shot> directory

        #defpending on the status of the run, some of these paths may not exist
        if cfg.clean >= 1:
            #top scripts
            flist = ["make_good_shots","make_hdrX.use","rback1","rback1_s","rback_field",
                    "rback_fix","rbacks","rbacks_2","rbfits","rbfits0","rbfits_fix","rbfits_s","rerun","rerun2","rfield",
                    "rfield.single","rfixspec","rimarb", "rsetdate","rsetdate0", "rtaras", "rtaremc",
                    "run0s", "run1s", "run2s", "runtar", "runtar0.slurm", "runtarm.defunct", "sun_use.dat",
                     "ifustat_20*.tab"]
            for fn in flist:
                system_command(cfg,f"rm {fn}")

            system_command(cfg,"unlink reductions")

            ############################
            #individual exposures
            # no do not remove) at this level ... we need these still if we need to rebuild the h5 file
            ############################
            #system_command(cfg,f"rm -rf d{cfg.datevshot.replace('v','s')}exp*")

            #########################
            #vdrp
            #########################
            if safe_cd(os.path.join(cfg.cwd,"vdrp")):
                #os.chdir(os.path.join(cfg.cwd,"vdrp"))
                #unlinks
                system_command(cfg, "unlink all.mch")
                system_command(cfg, "unlink norm.dat")
                system_command(cfg, "unlink shout.ifu")
                system_command(cfg, "unlink shout.ifustars")
                flist = glob.glob(f"{cfg.datevshot[0:8]}*")
                #exclude the .log
                flist.remove(f"{cfg.datevshot}.log")
                for fn in flist:
                    system_command(cfg,f"unlink {fn}")

                flist = glob.glob(f"radec*")
                for fn in flist:
                    system_command(cfg,f"unlink {fn}")

                system_command(cfg, "rm -rf match_pngs")

                system_command(cfg,"rm -rf fplane")
                system_command(cfg,"rm fplane*") #flat files
                system_command(cfg, "rm jan.txt")  # flat files
                system_command(cfg, "rm readme.install")  # flat files
                system_command(cfg, "rm setup.sh")  # flat files
                system_command(cfg, "rm astrometry.py")  # flat files

                system_command(cfg, "rm -rf vdrp")
                system_command(cfg, "rm vdrp*")
                system_command(cfg, "rm shuf*")

                if safe_cd(os.path.join(cfg.cwd, "vdrp/shifts")):
                    system_command(cfg, f"unlink dithall")
                    system_command(cfg, f"unlink {cfg.datevshot}")

                    flist = ["build_dithall","check_norms","check_shout.ifu","clean_dithall","clean_rta","cp2dithall",
                             "getpdf","getrdlist","gn.rt","jobsplitter","make_astrom_slurm","make_good_shots","mkhdr",
                             "plotxy.def","rbadsh","rbadsh0","rclean","rcleanshot","rgetdith","rhdr1.r1","rshotsrm","rta.*",
                             "runsh1","runsh2","run_shifts_old.sh","run_shifts.sh","tmp*"]
                    for fn in flist:
                        system_command(cfg,f"rm {fn}")




            ############################
            # alldet
            ############################
            if safe_cd(os.path.join(cfg.cwd,"alldet")):
                system_command(cfg, f"rm check*")
                system_command(cfg, f"rm make*")
                system_command(cfg, f"rm r*")
                system_command(cfg, f"rm fwhm.use")
                system_command(cfg, f"rm norm.dat")
                system_command(cfg, f"rm sky_res.use")

            ############################
            # cs
            ############################
            if safe_cd(os.path.join(cfg.cwd,"cs")):
                system_command(cfg, f"rm r*")
                system_command(cfg, f"rm make*")

            if safe_cd(os.path.join(cfg.cwd, "cs/spec")):
                system_command(cfg, f"rm r*")
                #system_command(cfg, f"rm *.tar") #moved to higher level clean
                #system_command(cfg, f"rm *.rcs")



            ####################################
            # extraneous spec from build h5?
            ####################################
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf spec")

            ##################################
            # detect
            ##################################
            if safe_cd(os.path.join(cfg.cwd,"detect")):
                system_command(cfg, f"rm -rf r_backup")
                system_command(cfg, f"rm -rf backup")
                system_command(cfg, f"rm -rf cal_script")
                system_command(cfg, f"rm -rf data")

                system_command(cfg, f"rm r*")
                system_command(cfg, f"rm make*")

                system_command(cfg,f"unlink {cfg.datevshot}.dithall")

                flist = ["update_fwhm_norm","check_tp","jobsplitter","extcor.dat","qfit_gaia.py","template.slurm","tpavg.dat"]
                for fn in flist:
                    system_command(cfg, f"rm {fn}")

                if safe_cd(os.path.join(cfg.cwd,f"detect/{cfg.datevshot}")):
                    flist = ["getsdssg","rfitfw0","rgetfw0", "rgettp0","rgettp0.karl", "rsp3fc","rspstar0","rspstar3",
                             "runstar0","j1","rgetfw","rgettp","rgettp0b","rsp3f","rspstar","rspstar2","runstar"]
                    for fn in flist:
                        system_command(cfg, f"rm {fn}")



            ####################################
            # Diagnose
            ####################################
            if safe_cd(cfg.cwd):
                system_command(cfg,"rm -rf ./Diagnose")
                flist = ["classification_000.fits","classification_001.fits","diagnose_spectra_cont.fits",
                         "diagnose_spectra_line.fits"]
                for fn in flist:
                    system_command(cfg,f"rm {fn}")



        if cfg.clean >= 2:  #more aggressive

            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm {cfg.datevshot.replace('v','')}_stats.pickle")
                system_command(cfg, f"rm {cfg.datevshot.replace('v','')}_ampstats.fits")

            #########################
            # vdrp
            #########################
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf vdrp")

            ##################################
            # detect
            ##################################
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf detect")
            #if safe_cd(os.path.join(cfg.cwd,"detect")):
            #    system_command(cfg, f"rm -rf {cfg.datevshot}")


            ############################
            # alldet
            ############################

            if safe_cd(os.path.join(cfg.cwd, "alldet")):
                system_command(cfg, f"rm -rf cal_out")
                system_command(cfg, f"rm -rf output")
                system_command(cfg, f"rm *.log")

            ############################
            # cs
            ############################
            if safe_cd(os.path.join(cfg.cwd,"cs")):
                system_command(cfg, f"rm *.log")

                if safe_cd(os.path.join(cfg.cwd, "cs/spec")):
                    system_command(cfg, f"rm *.log")
                    system_command(cfg, f"rm *.tar")
                    system_command(cfg, f"rm *.rcs")

            ##################################
            # ELiXer
            # can only clean this later as a negative clean value
            # since elixer is a manually queued slurm AFTER a run
            ##################################
            if cfg.clean_only and safe_cd(cfg.cwd):
                #elixer basic clean .. this leaves the top logs and the h5 files and all_pngs
                system_command(cfg, f"rm -rf ./elixer/line/dispatch_*")
                system_command(cfg, f"rm -rf ./elixer/cont/dispatch_*")
                system_command(cfg, f"rm ./elixer/line/ELIXER.o*")
                system_command(cfg, f"rm ./elixer/line/ELIXER.e*")
                system_command(cfg, f"rm ./elixer/line/elixer.run")
                system_command(cfg, f"rm ./elixer/line/elixer.slurm")
                system_command(cfg, f"rm ./elixer/cont/ELIXER.o*")
                system_command(cfg, f"rm ./elixer/cont/ELIXER.e*")
                system_command(cfg, f"rm ./elixer/cont/elixer.run")
                system_command(cfg, f"rm ./elixer/cont/elixer.slurm")


        if cfg.clean >= 3:  # still more aggressive
            #removes the last of the build analysis stuff and logs
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm {cfg.datevshot}_ampstats.tab")
                system_command(cfg, f"rm -rf analysis")
                system_command(cfg, f"rm -rf match_pngs")
                system_command(cfg, "rm *.log")

                ##################################
                # ELiXer
                # can only clean this later as a negative clean value
                # since elixer is a manually queued slurm AFTER a run
                ##################################
                if cfg.clean_only:
                    system_command(cfg, "rm elixer/*.log")
                    system_command(cfg, "rm elixer/*.dets")
                    system_command(cfg, "rm elixer/elixer_cmd.txt")
                    system_command(cfg, "rm elixer/line/*.log")
                    system_command(cfg, "rm elixer/cont/*.log")


        if cfg.clean >= 4:  # still more aggressive
            #removes the various original .mc, .spec, etc files
            #removed diagnose specific classifications
            # can no longer use --resume

            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf ./getcen")


            if safe_cd(cfg.cwd):
                ############################
                # individual exposures
                ############################
                system_command(cfg, f"rm -rf d{cfg.datevshot.replace('v', 's')}exp*")

                system_command(cfg, f"rm -rf alldet")
                system_command(cfg, f"rm -rf cs")



        if cfg.clean >= 5:  # ONLY keep the shot and detection h5 files
            if safe_cd(cfg.cwd):

                system_command(cfg, f"rm -rf elixer") #NOTE: here we can go ahead and remove since it is the top dir
                                                      #and does not matter if elixer has executed
                system_command(cfg, f"rm *.tab")
                system_command(cfg, f"rm *.pickle")
                system_command(cfg, f"rm *.fits")
                system_command(cfg, f"rm *.txt")
                system_command(cfg, f"rm *.dat")
                system_command(cfg, f"rm *.fwhm")
                system_command(cfg, f"rm status.*")

        cfg.clean_done = True

    except:
        print(f"Exception in post_clean(). {traceback.format_exc()}")


def progress_update(cfg,progress_dict,key=None,status=True):
    """
    update the progress record
    :param progress_dict:
    :return:
    """

    try:
        if key is not None:
            progress_dict[key] = status

        if cfg.write_progress:
            fn = os.path.join(cfg.cwd,"progress.dat")
            with open(fn, 'w') as f:
                json.dump(progress_dict, f, indent=4)  # indent=4 for pretty-printing

    except:
        print(f"Exception in progress_update(). {traceback.format_exc()}")

def progress_init(cfg,initialize_only=False):
    """
    read/initialize progress dict
    :param cfg:
    :return:
    """
    dtprog = None
    if not initialize_only:
        try:
            fn = os.path.join(cfg.cwd,"progress.dat")
            if os.path.exists(fn):
                with open(fn, 'r') as f:
                    dtprog = json.load(f)

                #assume this was successful, but also check --resume
                if cfg.resume is False:
                    print("********** NOTICE **********")
                    print("Call did NOT specfiy --resume, but this reduction appears to be at least partially complete.")
                    print("Will attempt to resume (implied) at the last incomplete step ...")
                    print("****************************")
                    cfg.resume = True

                if "s01b_amp_images" not in dtprog.keys():
                    dtprog["s01b_amp_images"] = False


                if "s04f_analysis" not in dtprog.keys():
                    dtprog["s04f_analysis"] = False

                    #consistency check / force (if ANY under a step are False, then the top of the step becomes False (incomplete))
                if "s04_make_shot" in dtprog.keys():
                    dtprog["s04_make_shot"] = dtprog["s04_make_shot"] and dtprog["s04a_get_ifucens"] and \
                                                    dtprog["s04b_rfft"] and dtprog["s04c_rcal_all"] and dtprog["s04d_shot_h5"] and \
                                                    dtprog["s04e_amp_stats"] and dtprog["s04f_analysis"]
                else:
                    dtprog["s04_make_shot"] = dtprog["s04_sky_subtraction"] and dtprog["s04a_get_ifucens"] and \
                                                    dtprog["s04b_rfft"] and dtprog["s04c_rcal_all"] and dtprog["s04d_shot_h5"] and \
                                                    dtprog["s04e_amp_stats"] and dtprog["s04f_analysis"]
                    dtprog.pop("s04_sky_subtraction") #remove the old name
                    #ordering is now wrong, but, programatically, it does not matter
                    #should be rare as we move forward, so just noting that the progress.dat for this will be
                    #out of order and move on



                dtprog["s05_detection"] = dtprog["s05_detection"] and dtprog["s05b_rdet_rf1"] and \
                                                dtprog["s05c_rgetmax"] and dtprog["s05e_detection_hdf5"]

                dtprog["s06_catalogs"] = dtprog["s06_catalogs"] and dtprog["s06b_fof"] and \
                                                dtprog["s06c_diagnose"] and dtprog["s06d_elixer"] and dtprog["s06e_source_cat"]
        except:
            print(f"Exception in progress_init(). {traceback.format_exc()}")

    if dtprog is None:
        dtprog = {"s01_run1s": False,
                  "s01b_amp_images" : False,
                  "s02_vdrp": False,
                  "s03_fluxcal": False,
                  "s04_make_shot": False,
                  "s04a_get_ifucens": False,
                  "s04b_rfft": False,
                  "s04c_rcal_all": False,
                  "s04d_shot_h5": False,
                  "s04e_amp_stats": False,
                  "s04f_analysis": False,
                  "s05_detection": False,
                  "s05b_rdet_rf1": False,
                  "s05c_rgetmax": False,
                  "s05e_detection_hdf5": False,
                  "s06_catalogs": False,
                  "s06b_fof": False,
                  "s06c_diagnose": False,
                  "s06d_elixer": False,
                  "s06e_source_cat": False, }

        if not initialize_only:
            progress_update(cfg,dtprog) #write out the file, no updates yet

        cfg.dtprog = dtprog
    return dtprog



def write_summary(cfg):
    """

    :param cfg:
    :return:
    """
    orig_dir = os.getcwd()
    try:
        os.chdir(cfg.cwd)
        if not os.path.exists(f"{cfg.datevshot}.h5"):
            write_limited_summary(cfg)
            os.chdir(orig_dir)
            return

        h5 = tables.open_file(f"{cfg.datevshot}.h5",mode="r")

        cfg.set_warn = False

        try:
            ra = h5.root.Shot.read(field='ra')[0]
            dec = h5.root.Shot.read(field='dec')[0]
        except:
            ra = -999.9
            dec = -999.9

        with open("summary.txt","w") as f:
            f.write(f"shot:\t\t{cfg.datevshot}\n")
            f.write(f"field:\t\t{h5.root.Shot.read(field='field')[0]} - {cfg.shot_obj}\n")
            f.write(f"RA,Dec:\t\t{ra}, {dec}\n")
            f.write(f"exptimes:\t{h5.root.Shot.read(field='exptime')[0]}\n")
            dit_norms = h5.root.Shot.read(field="relflux_virus")[0]
            f.write(f"relflux_virus:\t{dit_norms}\n")
            try:
                max_ditnorm = np.max(dit_norms) / np.min(dit_norms)
                f.write(f"dither norm:\t{max_ditnorm:0.4f}\n")
            except:
                cfg.set_warn = True

            try:
                waves = h5.root.Calibration.Throughput.throughput.read(field="wavelength")
                tp = h5.root.Calibration.Throughput.throughput.read(field="throughput")
                tp_at_w = np.interp(4600.0, waves, tp)
            except:
                tp_at_w = -999.9
                cfg.set_warn = True

            try:
                if ra > -999:
                    coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
                    dust_corr = deredden_spectra([4540.], coord)[0]
                    f.write(f"dust_corr@4540:\t{dust_corr:0.2f}x\n")
                else:
                    f.write(f"dust_corr@4540:\t???\n")
            except:
                print(f"[{cfg}] Error computing dust correction for shot.", traceback.format_exc())
                f.write(f"dust_corr@4540:\t???\n")

            try:
                f.write(f"response_4540:\t{h5.root.Shot.read(field='response_4540')[0]:0.4f}\n")
            except:
                f.write(f"response_4540:\t???\n")
                cfg.set_warn = True

            try:
                f.write(f"interp  @4600:\t{tp_at_w:0.4f}\n")
            except:
                f.write(f"interp  @4600:\t???\n")
                cfg.set_warn = True

            try:
                f.write(f"fwhm_virus:\t{h5.root.Shot.read(field='fwhm_virus')[0]:0.4f}\n")
            except:
                f.write(f"fwhm_virus:\t???\n")
                cfg.set_warn = True

            try:
                dx1 = h5.root.Shot.read(field="xditherpos")[0]
                dy1 = h5.root.Shot.read(field="yditherpos")[0]
                f.write(f"Dithers: \t({dx1[0]:0.4f},{dy1[0]:0.4f}) ({dx1[1]:0.4f},{dy1[1]:0.4f}) ({dx1[2]:0.4f},{dy1[2]:0.4f})\n")
            except:
                f.write(f"Dithers: \t???\n")
                cfg.set_warn = True


            if cfg.avg_sky is None and cfg.resume: #this is true on a resume
                cfg.avg_sky = get_avg_sky(cfg)

            if cfg.avg_sky is not None:
                f.write(f"Avg sky: \t{cfg.avg_sky}\n")
            else:
                f.write(f"Avg sky: \t???\n")
                cfg.set_warn = True

            try:
                total = 0
                with open(f"{cfg.cwd}/elixer/line.dets", "r") as f2:
                    total = len(f2.readlines())

                #outer file:
                f.write(f"Line Dets:\t{total}\n")

            except:
                f.write(f"Line Dets:\t???\n")

            try:
                total = 0
                with open(f"{cfg.cwd}/elixer/cont.dets", "r") as f2:
                    total = len(f2.readlines())

                # outer file:
                f.write(f"Cont Dets:\t{total}\n")

            except:
                f.write(f"Cont Dets:\t???\n")

            try:
                # outer file:

                if cfg.num_bad_amps < 0:
                    _ = count_amps(cfg)

                f.write(f"Bad Amps:\t{cfg.num_bad_amps}/{cfg.num_all_amps}\n")
                if cfg.amp_stats_problem > 0:
                    f.write(f"Amp Stats:\tFail ({cfg.amp_stats_problem})\n")
                    cfg.set_warn = True
                elif cfg.amp_stats_problem < 0:
                    f.write(f"Amp Stats:\tUnspecified\n")
                else:
                    f.write(f"Amp Stats:\tPass\n")
            except:
                f.write(f"Bad Amps:\t???\n")
                f.write(f"Amp Stats:\t???\n")

            f.write("") #always end with a blank line, so can read easier with cat
    except:
        print(f"[{cfg.datevshot}] Exception! trying to write summary file", traceback.format_exc())

    try:
        h5.close()
    except:
        pass

    os.chdir(orig_dir)


def write_limited_summary(cfg):
    """
    like write summary, but with whatever we have
    normally if this fails to get to the .h5 file creation

    :param cfg:
    :return:
    """
    orig_dir = os.getcwd()
    try:
        os.chdir(cfg.cwd)

        cfg.set_warn = False

        #pre-check ... the summary might already exist and be populated
        #  and if this is an aborted or limited run, do not overwrite that
        #  we'll use reflux_virus as a good indicator that the summary has real data in it
        #    and data that this limited write would not be able to update
        overwrite = True
        if os.path.exists("summary.txt"):
            overwrite = False
            with open("summary.txt","r") as f:
                for line in f:
                    if "relflux_virus" in line and "[]" in line:
                        #this is not populated
                        overwrite = True
                        break

        if not overwrite:
            print(f"[{cfg.datevshot}] existing, populated summary.txt file found. Will not overwrite.")
            return

        if cfg.shot_ra == -999.9:
            get_local_ra_dec(cfg)

        try:
            ra = cfg.shot_ra
            dec = cfg.shot_dec
        except:
            ra = -999.9
            dec = -999.9

        # if ra == -999.9: #there is probably nothing available to write
        #     os.chdir(orig_dir)
        #     return

        with open("summary.txt","w") as f:
            f.write(f"shot:\t\t{cfg.datevshot}\n")
            f.write(f"field:\t\t{cfg.shot_obj}\n")
            f.write(f"RA,Dec:\t\t{ra:0.6f}, {dec:0.6f}\n")
            f.write(f"exptimes:\t{cfg.total_exp_time:0.1f}s total\n")
            dit_norms = cfg.dither_norms
            f.write(f"relflux_virus:\t{dit_norms}\n")
            try:
                max_ditnorm = np.max(dit_norms) / np.min(dit_norms)
                f.write(f"dither norm:\t{max_ditnorm:0.4f}\n")
            except:
                cfg.set_warn = True

            # try:
            #     waves = h5.root.Calibration.Throughput.throughput.read(field="wavelength")
            #     tp = h5.root.Calibration.Throughput.throughput.read(field="throughput")
            #     tp_at_w = np.interp(4600.0, waves, tp)
            # except:
            #     tp_at_w = -999.9
            #     cfg.set_warn = True

            try:
                if ra > -999:
                    coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
                    dust_corr = deredden_spectra([4540.], coord)[0]
                    f.write(f"dust_corr@4540:\t{dust_corr:0.2f}x\n")
                else:
                    f.write(f"dust_corr@4540:\t???\n")
            except:
                print(f"[{cfg}] Error computing dust correction for shot.", traceback.format_exc())
                f.write(f"dust_corr@4540:\t???\n")

            try:
                f.write(f"response_4540:\t???\n")
            except:
                f.write(f"response_4540:\t???\n")
                cfg.set_warn = True

            try:
                f.write(f"interp  @4600:\t???\n")
            except:
                f.write(f"interp  @4600:\t???\n")
                cfg.set_warn = True

            try:
                f.write(f"fwhm_virus:\t???\n")
            except:
                f.write(f"fwhm_virus:\t???\n")
                cfg.set_warn = True

            try:
                f.write(f"Dithers: \t({cfg.dither_positions})\n")
            except:
                f.write(f"Dithers: \t???\n")
                cfg.set_warn = True


            if cfg.avg_sky is None and cfg.resume: #this is true on a resume
                cfg.avg_sky = get_avg_sky(cfg)

            if cfg.avg_sky is not None:
                f.write(f"Avg sky: \t{cfg.avg_sky}\n")
            else:
                f.write(f"Avg sky: \t???\n")
                cfg.set_warn = True

            try:
                total = 0
                with open(f"{cfg.cwd}/elixer/line.dets", "r") as f2:
                    total = len(f2.readlines())

                #outer file:
                f.write(f"Line Dets:\t{total}\n")

            except:
                f.write(f"Line Dets:\t???\n")

            try:
                total = 0
                with open(f"{cfg.cwd}/elixer/cont.dets", "r") as f2:
                    total = len(f2.readlines())

                # outer file:
                f.write(f"Cont Dets:\t{total}\n")

            except:
                f.write(f"Cont Dets:\t???\n")

            try:
                # outer file:

                if cfg.num_bad_amps < 0:
                    _ = count_amps(cfg)

                f.write(f"Bad Amps:\t{cfg.num_bad_amps}/{cfg.num_all_amps}\n")
                if cfg.amp_stats_problem > 0:
                    f.write(f"Amp Stats:\tFail ({cfg.amp_stats_problem})\n")
                    cfg.set_warn = True
                elif cfg.amp_stats_problem < 0:
                    f.write(f"Amp Stats:\tUnspecified\n")
                else:
                    f.write(f"Amp Stats:\tPass\n")
            except:
                f.write(f"Bad Amps:\t???\n")
                f.write(f"Amp Stats:\t???\n")

            f.write("")  # always end with a blank line, so can read easier with cat
    except:
        print(f"[{cfg.datevshot}] Exception! trying to write limited summary file", traceback.format_exc())


    os.chdir(orig_dir)


def Quit(cfg,rc,msg=None,do_write_status=True,do_post_clean=True,do_write_summary=True):
    """

    :param cfg:
    :param rc:
    :param msg:
    :return:
    """

    if msg is not None:
        print(f"[{cfg.datevshot}] ({rc})",msg)
    else:
        print(f"[{cfg.datevshot}] ({rc})")

    dtprog = cfg.dtprog
    if dtprog is None:
        dtprog = progress_init(cfg,initialize_only=True)
        cfg.write_progress = False

    #on a resume, make sure we have loaded the warn triggers
    if cfg.resume and (do_write_status or do_write_summary):
        if cfg.amp_stats_problem < 0 and dtprog["s04e_amp_stats"]:
            #this is a resume, we have not set amp_stats_problem AND at some point in the last run, amp_stats were computed
            #so we probably skipped it this time and need to know this information for below
            #this will set amp_stats_problem, num_bad_amps, bnum_all_amps, etc
            amp_stats(cfg,update=False)

        #rfft needs to have already run
        if cfg.avg_sky is None and dtprog["s04b_rfft"]:
            avg_sky = get_avg_sky(cfg)
            if avg_sky is not None:
                cfg.avg_sky = avg_sky
                if avg_sky > MAX_SAFE_AVG_SKY:
                    if avg_sky >= FAIL_AVG_SKY:
                        print(
                            f"[{cfg.datevshot}] Average Sky is catastrophically large and/or unable to be fit: {avg_sky:0.1f}.")
                        rc = -1  # fatal
                    else:
                        # not fatal, but still a warning
                        print(f"[{cfg.datevshot}] Average Sky is problematically large: {avg_sky:0.1f}. Non-fatal.")

        if (cfg.num_cont_dets < 0 or cfg.num_line_dets < 0) and dtprog["s06b_fof"]:
            get_detection_counts(cfg)


    #?always write the summary??
    if do_write_summary:
        write_summary(cfg)

    #if rc >=0 and not cfg.multifits_only:
    #    write_summary(cfg)

   # if cfg.orig_stdout:
   #     sys.stdout = cfg.orig_stdout

    #if cfg.orig_stderr:
    #    sys.stderr = cfg.orig_stderr

    if cfg.file_stdout:
        cfg.file_stdout.flush()
        cfg.file_stdout.close()

        #repeat the final message to the console
       # if msg is not None:
        #    print(f"({rc})", msg)
        #else:
        #    print(f"({rc})")

    try: #remove the status.run file (if it exists)
        if os.path.exists(os.path.join(cfg.cwd,'status.run')):
            system_command(cfg,f"rm {os.path.join(cfg.cwd,'status.run')}")
    except:
        pass

    try:
        if do_write_status and safe_cd(cfg.cwd):
            status_str = "status"

            #
            # note: each status line should be of the form
            #   [<datevsshot>] <status>. <msg>

            if os.path.exists(f"{status_str}.fail") or \
               os.path.exists(f"{status_str}.warn") or \
               os.path.exists(f"{status_str}.pass"):
                #this is probably a resume
                if cfg.resume:
                    status_str = "status.resume"
                    #we will report the status of just this bit, but leave the original one around
                    #most likely you have a "warn" that was partly re-run and that re-run bit was fine
                    #  so you'd get a status.warn and a status.resume.pass
                else:
                    #wipe out the old one and replace
                    system_command(cfg, f"rm {status_str}.*")

            if rc < 0:
                with open(f"{status_str}.fail","w") as f:
                    f.write(f"[{cfg.datevshot}] fail. ({rc}) {msg}\n")
            else:

                if cfg.num_line_dets > WARN_NUM_LINE_DETS:
                    cfg.set_warn = True
                if cfg.num_cont_dets > WARN_NUM_CONT_DETS:
                    cfg.set_warn = True
                if cfg.avg_sky is not None and cfg.avg_sky > MAX_SAFE_AVG_SKY:
                    cfg.set_warn = True

                if cfg.num_bad_amps < 0:
                    _ = count_amps(cfg)

                if cfg.num_bad_amps > 0 and cfg.num_all_amps > 0:
                    if cfg.num_bad_amps / cfg.num_all_amps >= 0.5:
                        #half or more are bad
                        cfg.set_warn = True
                        if cfg.amp_stats_problem <= 0:
                            cfg.amp_stats_problem = 1


                if cfg.set_warn:
                    #note: if avg_sky > FAIL_AVG_SKY, then the failure condition would have already tripped
                    #and we would be in the case above (rc < 0)
                    with open(f"{status_str}.warn", "w") as f:
                        f.write(f"[{cfg.datevshot}] warn. ({rc}) {msg}\n")
                        if (cfg.avg_sky is not None and cfg.avg_sky > MAX_SAFE_AVG_SKY):
                            f.write(f"[{cfg.datevshot}] warn. Avg Sky is large: {cfg.avg_sky:0.1f} \n")
                        if cfg.amp_stats_problem != 0:
                            f.write(f"[{cfg.datevshot}] warn. Amp Stats Issue ({cfg.amp_stats_problem}),"
                                    f" {cfg.num_bad_amps}/{cfg.num_all_amps} bad amps\n")
                        if cfg.num_line_dets > WARN_NUM_LINE_DETS:  # pretty generous here, maybe should also adjust with exptime and conditions?
                            f.write(f"[{cfg.datevshot}] warn. Num of line dets is large: {cfg.num_line_dets} \n")
                        if cfg.num_cont_dets > WARN_NUM_CONT_DETS:
                            f.write(f"[{cfg.datevshot}] warn. Num of cont dets is large: {cfg.num_cont_dets} \n")
                else:
                    with open(f"{status_str}.pass","w") as f:
                        f.write(f"[{cfg.datevshot}] pass. ({rc}) {msg}\n")
    except:
        pass

    node_clean(cfg) #always clean the node

    if rc >= 0 and do_post_clean:
        post_clean(cfg)

    exit(rc)


def blocking_command(cfg,cmd):
    """
    wrapper to subprocess Popen
    blocks and returns result

    :param cfg:
    :param cmd:
    :return:
    """

    try:
        rc = 0
        #cmdlist = [cmd[0], cmd[1:]]
        if EchoCmds:
            print(f"[{cfg.datevshot}] subproc CMD: ({os.getcwd()}) > {cmd}")

        rc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

        if EchoCmds:
            print(f"[{cfg.datevshot}] subproc rc = {rc}: CMD = ({os.getcwd()}) > {cmd}")
    except:
        print(f"[{cfg.datevshot}] Exception! in blocking_command. cmd = {cmd}\n",traceback.format_exc())
        rc = -1

    return rc

def system_command(cfg,cmd):
    """

    wrapper to execute a system command

    :param cfg:
    :param cmd:
    :return:
    """

    #echo the command
    if EchoCmds:
        print(f"[{cfg.datevshot}] CMD: ({os.getcwd()}) > {cmd}")

    if cfg.file_stdout:
        os.system(f"{cmd} &>> {cfg.file_stdout.name}")
    else:
        os.system(f"{cmd}")



def linear_exptime_scale(cfg):
    """
    based on HETDEX normalization on exposure time, return a linear scaling multiplier

    :param cfg:
    :return:
    """
    mux = 1.0
    try:
        if cfg.total_exp_time is not None and cfg.numexp is not None:
            mux = max(1.0, (cfg.total_exp_time / max(1, cfg.numexp) / 360.))  # 360secs is the nominal default for a HETDEX exposure
    except:
        print(f"[{cfg.datevshot}] linear_exptime_scale excetpion. {traceback.format_exc()}")

    return mux

def get_exposure_times(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        date = cfg.datevshot[0:8]
        path = os.path.join(HET_by_date, date)
        virus_shot = "virus0000" + cfg.datevshot[-3:]
        #base_tarfn = os.path.join(path, f"virus/{virus_shot}.tar")
        if cfg.virus_tar_path is not None:
            base_tarfn = cfg.virus_tar_path
            if base_tarfn[-4:] != ".tar": #this is probably an intermediate directory
                extended_tarfn = os.path.join(base_tarfn,f"{cfg.datevshot[0:8]}/virus/{virus_shot}.tar")
                if os.path.exists(extended_tarfn):
                    base_tarfn = extended_tarfn
                else:
                    extended_tarfn = os.path.join(base_tarfn, f"/virus/{virus_shot}.tar")
                    if os.path.exists(extended_tarfn):
                        base_tarfn = extended_tarfn
                    else:
                        extended_tarfn = os.path.join(base_tarfn, f"{virus_shot}.tar")
                        if os.path.exists(extended_tarfn):
                            base_tarfn = extended_tarfn
                        else: #did not find it ... could try to scan for any {virus_shot}.tar below here
                            print(f"[{cfg.datevshot}] could not find tarfile?")
        else:
            base_tarfn = os.path.join(path, f"virus/{virus_shot}.tar")

        exposure_times = []
        if os.path.exists(base_tarfn):
            with tar.open(base_tarfn, "r") as tarfh:
                fns = np.array(tarfh.getnames())
                #only want the files that are virus*/exp??/virus/....
                fns_sel = np.array(["/exp" in fn for fn in fns])
                fns = fns[fns_sel]
                fns_sel = np.array([".fits" in fn for fn in fns])
                fns = fns[fns_sel]
                # should look like a list of:  virus0000007/exp01/virus/20241017T024257.2_106RU_sci.fits
                # should all be the same base within each exposure
                # not all the file names are that way, though

                all_exps = np.array([x.split("/")[1] for x in fns])
                exps = np.unique(all_exps)
                # just need one from each exp
                for exp in exps:
                    sel = all_exps == exp
                    if np.count_nonzero(sel) > 0:
                        fn0 = fns[sel][0]
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, fn0)
                        fh = fits.open(t1)

                        #there is an apparent issue in populating the EXPTIME for some observations
                        # perhaps related to the use of the HPF
                        try:
                            if abs(fh[0].header['EXPTIME'] - fh[0].header['PEXPTIME']) > 30.0:
                                use_exp = fh[0].header['PEXPTIME'] - 8.0 #the 8.0 is an approximate correction from Greg Z.
                                print(f"[{cfg.datevshot}] EXPTIME {fh[0].header['EXPTIME']} vs "
                                      f"PEXPTIME {fh[0].header['PEXPTIME']} large difference. Using PEXPTIME - 8.0s.")
                            else:
                                use_exp = fh[0].header['EXPTIME']
                        except:
                            use_exp = fh[0].header['EXPTIME']
                            print(f"[{cfg.datevshot}] Exception retrieving EXPTIME. Using value: {use_exp}",
                                  traceback.format_exc())

                        exposure_times.append(use_exp)
                        fh.close()

        if len(exposure_times) == 0:
            # could not find any
            print(f"[{cfg.datevshot}] Total exposure time: fail, could not compute")
            cfg.total_exp_time = 0.0
        else:
            cfg.total_exp_time = np.nansum(exposure_times)
            cfg.dither_exp_times = exposure_times
            print(f"[{cfg.datevshot}] Total exposure time: {cfg.total_exp_time}")

    except:
        print(f"[{cfg.datevshot}] Exception in get_exposure_times: {traceback.format_exc()}")

def get_guider_fwhm(cfg):
    """
    try to get the seeing fwhm (IQ -- image quality) from the guider (gc1 or gc2) for the shot
    :param cfg:
    :return:
    """


    def get_gc_path(cfg,which_gc="gc1",path=""):
        base_tarfn = os.path.join(path,f"{which_gc}/{which_gc}.tar")
        if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
            #try based on the virus_tar_path
            base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"),f"{which_gc}.tar")

            if not os.path.exists(base_tarfn):
                # try based on the virus_tar_path, up two levels
                base_tarfn = os.path.join(os.path.dirname(cfg.virus_tar_path), f"../{which_gc}/{which_gc}.tar")

                if not os.path.exists(base_tarfn):
                    #sometimes is called images.tar
                    base_tarfn = os.path.join(os.path.dirname(cfg.virus_tar_path), f"../{which_gc}/images.tar")

                    if not os.path.exists(base_tarfn):
                        base_tarfn = f"./{which_gc}.tar" #last chance, check locally

        return base_tarfn

    def which_dither_time(curr_time="",dither_starts=[]):
        #starts look like yyyymmddThhmmss.s
        #can use alphabet sorting
        try:
            idx = sorted(list(dither_starts) + [curr_time]).index(curr_time)
        except:
            idx = 99
        return idx if idx > 0 else 1

    try:
        saved_fn = "guider.fwhm"
        saved_detailed_fn = "detailed.fwhm" #each entry from the guider over the whole shot, all dithers, all times
                                            #this is extra info for later, optional diagnostics
                                            #and is only gathered if GuiderFWHM_ALL is set

        if os.path.exists(saved_fn):
            fwhm = np.loadtxt(saved_fn,dtype=float) #just one value
            if fwhm is None or np.isnan(fwhm):
                 #we will continue below and rebuild
                print(f"[{cfg.datevshot}] Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
            else:
                fail = False
                try:
                    if len(fwhm) == 0 or fwhm[0] <= 0:
                        print(f"[{cfg.datevshot}] Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
                        fail = True
                except:
                    pass

                try:
                    fwhm = float(fwhm)
                except:
                    print(f"[{cfg.datevshot}] Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
                    fail = True

                if not fail:
                    print(f"[{cfg.datevshot}] Found {saved_fn}. Using guider seeing FWHM = {fwhm}")
                    return fwhm
        else:
            print(f"[{cfg.datevshot}] Did not find {saved_fn}.")

        exposure_times = [] #exposure times (seconds) from HDU
        exposure_fn_times =[] #date T time string from the fileanames
        #dither_start_times = []  # same time as exposure_fn_times  ... the time on the file name matches the UT time the exposure STARTED

        gc_start_times = []
        gc_stop_times = []


        gc1_names = None
        gc2_names = None
        #gc1_times = None
        #gc2_times = None

        gc1_near = None
        gc2_near = None

        date = cfg.datevshot[0:8]
        path = os.path.join(HET_by_date,date)

        #first what is the virus shot(s) we need
        #the guider filenames will not be the same, but will be close
        virus_shot = "virus0000" + cfg.datevshot[-3:]

        if cfg.virus_tar_path is not None:
            base_tarfn = cfg.virus_tar_path
            if os.path.isdir(base_tarfn):
                base_tarfn = os.path.join(cfg.virus_tar_path,f"{cfg.datevshot[0:8]}/virus/{virus_shot}.tar")
        else:
            base_tarfn = os.path.join(path, f"virus/{virus_shot}.tar")

        if os.path.exists(base_tarfn):
            with tar.open(base_tarfn, "r") as tarfh:
                fns = np.array(tarfh.getnames())
                #should look like a list of:  virus0000007/exp01/virus/20241017T024257.2_106RU_sci.fits
                #should all be the same base within each exposure
                #only want the files that are virus*/exp??/virus/....
                fns_sel = np.array(["/exp" in fn for fn in fns])
                fns = fns[fns_sel]
                fns_sel = np.array([".fits" in fn for fn in fns])
                fns = fns[fns_sel]

                all_exps = np.array([x.split("/")[1] for x in fns])
                exps = np.unique(all_exps)
                #just need one from each exp
                for exp in exps:
                    sel = all_exps == exp
                    if np.count_nonzero(sel) > 0:
                        fn0 = fns[sel][0]
                        name = fn0.split("/")[-1].split("_")[0]
                        exposure_fn_times.append(name)

                        t1, p1 = Utils.open_file_from_tar(base_tarfn, fn0)
                        fh = fits.open(t1)

                        #there is an apparent issue in populating the EXPTIME for some observations
                        # perhaps related to the use of the HPF
                        try:
                            if abs(fh[0].header['EXPTIME'] - fh[0].header['PEXPTIME']) > 30.0:
                                use_exp = fh[0].header['PEXPTIME'] - 8.0  # the 8.0s is an approximate correction from Greg Z.
                                print(f"[{cfg.datevshot}] EXPTIME {fh[0].header['EXPTIME']} vs "
                                      f"PEXPTIME {fh[0].header['PEXPTIME']} large difference. Using PEXPTIME - 8.0s.")
                            else:
                                use_exp = fh[0].header['EXPTIME']
                        except:
                            use_exp = fh[0].header['EXPTIME']
                            print(f"[{cfg.datevshot}] Exception retrieving EXPTIME. Using value: {use_exp}",
                                  traceback.format_exc())

                        exposure_times.append(use_exp)
                        #exposure_times.append(fh[0].header['EXPTIME'])
                        #dither_start_times.append(fh[0].header['UT'])  #note: 'DATE' is when the file was written out
                        #or would DARKTIME be better ??
                        fh.close()


        if len(exposure_fn_times) == 0:
            #could not find any
            print(f"[{cfg.datevshot}] Total exposure time: fail, could not locate")
            return None
        else:
            cfg.total_exp_time = np.nansum(exposure_times)
            cfg.dither_exp_times = exposure_times
            print(f"[{cfg.datevshot}] Total exposure time: {cfg.total_exp_time}")

        #try gc1
        # base_tarfn = os.path.join(path,"gc1/gc1.tar")
        # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
        #     #try based on the virus_tar_path, up two levels
        #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"),"gc1.tar")

        print(f"[{cfg.datevshot}] Computing guider for seeing FWHM (this can take a while) ...")

        base_tarfn = get_gc_path(cfg,"gc1",path)
        if os.path.exists(base_tarfn):
            print(f"[{cfg.datevshot}] Getting gc1 from {base_tarfn}")
            with tar.open(base_tarfn, "r") as tarfh:
                gc1_names = tarfh.getnames()
                #should look like a list of:  20241017T024256.0_gc1_sci.fits
                #only want the dateTtime prefix
                #gc1_times = [x.split("/")[1].split("_")[0] for x in gc1_names]



        #try gc2
        # base_tarfn = os.path.join(path,"gc2/gc2.tar")
        # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
        #
        #     #try top most:
        #
        #
        #     #then try under the virusXXXXXX, but these don't have the IQ card?
        #     #try based on the virus_tar_path, up two levels
        #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"),"gc2.tar")

        base_tarfn = get_gc_path(cfg, "gc2", path)
        if os.path.exists(base_tarfn):
            print(f"[{cfg.datevshot}] Getting gc2 from {base_tarfn}")
            with tar.open(base_tarfn, "r") as tarfh:
                gc2_names = tarfh.getnames()
                #should look like a list of:  20241017T024256.0_gc1_sci.fits
                #only want the dateTtime prefix ... actually, can just use as is
                #gc2_times = [x.split("/")[1].split("_")[0] for x in gc2_names]


        # maybe get the exposure time for the science image and then get all the
        # guider images that overlap with that time?
        # get the time from the filename, subtract the exposure time and accept all
        # guider images that have a filename time AFTER and within some small XX beyond the original sciece filetime


        if GuiderFWHM_ALL: #use all within specified time range

            print(f"[{cfg.datevshot}] Getting approx start/stop times for science frames ...")
            for fntime, exptime in zip(exposure_fn_times, exposure_times):
                #fntime, from the name is the time the exposure STARTED (not written). It matches to "UT" card and the stop is the "DATE" card in HDU
                fntime = datetime(int(fntime[0:4]), int(fntime[4:6]), int(fntime[6:8]),
                                  int(fntime[9:11]), int(fntime[11:13]), int(fntime[13:15]))
                delta_time = fntime #- timedelta(seconds=6.0) # !!! YES this is correct. The filename has the START time, not the writeout (stop)time
                #delta_time = fntime - timedelta(seconds=exptime)
                time_str = f"{str(delta_time.year)}{str(delta_time.month).zfill(2)}{str(delta_time.day).zfill(2)}T" \
                           f"{str(delta_time.hour).zfill(2)}{str(delta_time.minute).zfill(2)}{str(delta_time.second).zfill(2)}.0"

                gc_start_times.append(time_str)

                #Guider typically records every few seconds, so no need to go after
                #could go an extra one past just to get one more data point 5 -6 seconds is typical
                delta_time = fntime + timedelta(seconds=exptime) #+ timedelta(seconds=6.0)   #timedelta(seconds=120.0)  # go up to 2 minutes after
                #delta_time = fntime
                time_str = f"{str(delta_time.year)}{str(delta_time.month).zfill(2)}{str(delta_time.day).zfill(2)}T" \
                           f"{str(delta_time.hour).zfill(2)}{str(delta_time.minute).zfill(2)}{str(delta_time.second).zfill(2)}.0"
                gc_stop_times.append(time_str)

            if gc1_names is not None:
                gc1_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                print(f"[{cfg.datevshot}] Checking corresponding {len(gc1_names)} gc1 files and {len(gc_start_times)} "
                      f"timestamps ...")
                for start_time, stop_time in zip(gc_start_times,gc_stop_times):
                    for name in gc1_names:
                        #just want the time part
                        try:
                            time_part = os.path.basename(name)[:17] #[2:19]
                            if start_time <= time_part <= stop_time: #yes, technically string compares, but works for the time format here
                                gc1_near.append(name)
                        except:
                            pass

            if gc2_names is not None:
                gc2_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                print(f"[{cfg.datevshot}] Checking corresponding {len(gc2_names)} gc2 files and {len(gc_start_times)} "
                      f"timestamps ...")
                for start_time, stop_time in zip(gc_start_times,gc_stop_times):
                    for name in gc2_names:
                        #just want the time part
                        try:
                            time_part = os.path.basename(name)[:17]#[2:19]
                            if start_time <= time_part <= stop_time: #yes, technically string compares, but works for the time format here
                                gc2_near.append(name)
                        except:
                            pass
        else: #just use the two nearest
            # get nearest (sort with exposusre name, find the index of the exposure name and take the index before and after?)

            if gc1_names is not None:
                gc1_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                for expname in exposure_fn_times:
                    fake_name = f"./{expname}_gc1_sci.fits"
                    gc_x = sorted(gc1_names + [fake_name])
                    idx = gc_x.index(fake_name)
                    lidx = max(0,idx-1)
                    ridx = min(len(gc_x),idx+1)
                    #could be an exact match, but if so, either left or right will then also be the same
                    gc1_near.append([gc_x[lidx]] + [gc_x[ridx]])

            if gc2_names is not None:
                gc2_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                for expname in exposure_fn_times:
                    fake_name = f"./{expname}_gc2_sci.fits"
                    gc_x = sorted(gc2_names + [fake_name])
                    idx = gc_x.index(fake_name)
                    lidx = max(0,idx-1)
                    ridx = min(len(gc_x),idx+1)
                    #could be an exact match, but if so, either left or right will then also be the same
                    gc2_near.append([gc_x[lidx]] + [gc_x[ridx]])

        iq = [] #image quality (seeing fwhm) list
        iq_time = [] #time stamp for the image quality
        iq_dither = []

        if GuiderFWHM_ALL:
            #note: which guider is active could change? so check both?
            # base_tarfn = os.path.join(path, "gc1/gc1.tar")
            # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
            #     # try based on the virus_tar_path, up two levels
            #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc1.tar")

            base_tarfn = get_gc_path(cfg, "gc1", path)

            excessive_time_count = 0
            iq_idx_start =0
            if os.path.exists(base_tarfn) and gc1_near is not None and len(gc1_near) > 0:
                print(f"[{cfg.datevshot}] Collecting data from matched {len(gc1_near)} gc1 files ...")
                for name in gc1_near:
                #for name in tqdm(gc1_near):
                    try:
                        start_time = time.perf_counter_ns()

                        t1, p1 = Utils.open_file_from_tar(base_tarfn, name)
                        fh = fits.open(t1)
                        if fh[0].header['GUIDLOOP'] == 'ACTIVE':
                            try:
                                iq.append(float(fh[0].header['IQ'])) #this is clipped later to exlcude 0s and other bad values
                                iq_time.append(os.path.basename(name)[:17])
                                iq_dither.append(which_dither_time(iq_time[-1],exposure_fn_times))
                            except:
                                try:
                                    #card not there, try using PIXSCALE and OBFWHM
                                    pixscale = float(fh[0].header['PIXSCALE'])
                                    obfwhm = float(fh[0].header['OBFWHM'])
                                    iq.append(obfwhm * pixscale)
                                    iq_time.append(os.path.basename(name)[:17])
                                    iq_dither.append(which_dither_time(iq_time[-1],exposure_fn_times))
                                except:
                                    pass

                        fh.close()
                        t1.close()
                        elapsed = (time.perf_counter_ns() - start_time) // 1e9
                        if elapsed > 5.0: #5 seconds is excessive
                            if excessive_time_count > 2:
                                print(f"[{cfg.datevshot}] Excessive time collecting gc1 data. Aborting collection.")
                                iq = iq[0:iq_idx_start + 1]
                                iq_time = iq_time[0:iq_idx_start + 1]
                                iq_dither = iq_dither[0:iq_idx_start + 1]
                                break
                            else:
                                excessive_time_count += 1
                        else:
                            excessive_time_count = 0


                    except:
                        pass #don't let one bad read bomb out

            # base_tarfn = os.path.join(path, "gc2/gc2.tar")
            # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
            #     # try based on the virus_tar_path, up two levels
            #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc2.tar")

            base_tarfn = get_gc_path(cfg, "gc2", path)
            excessive_time_count = 0
            iq_idx_start = len(iq)
            if os.path.exists(base_tarfn) and gc2_near is not None and len(gc2_near) > 0:
                print(f"[{cfg.datevshot}] Collecting data from matched {len(gc2_near)} gc2 files ...")
                for name in gc2_near:
                #for name in tqdm(gc2_near):
                    try:
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, name)
                        fh = fits.open(t1)
                        if fh[0].header['GUIDLOOP'] == 'ACTIVE':
                            try:
                                iq.append(float(fh[0].header['IQ'])) #this is clipped later to exlcude 0s and other bad values
                                iq_time.append(os.path.basename(name)[:17])
                                iq_dither.append(which_dither_time(iq_time[-1],exposure_fn_times))
                            except:
                                try:
                                    # card not there, try using PIXSCALE and OBFWHM
                                    pixscale = float(fh[0].header['PIXSCALE'])
                                    obfwhm = float(fh[0].header['OBFWHM'])
                                    iq.append(obfwhm * pixscale)
                                    iq_time.append(os.path.basename(name)[:17])
                                    iq_dither.append(which_dither_time(iq_time[-1],exposure_fn_times))
                                except:
                                    pass

                        fh.close()
                        t1.close()
                        elapsed = (time.perf_counter_ns() - start_time) // 1e9
                        if elapsed > 5.0: #5 seconds is excessive
                            if excessive_time_count > 2:
                                print(f"[{cfg.datevshot}] Excessive time collecting gc2 data. Aborting collection.")
                                iq = iq[0:iq_idx_start + 1]
                                iq_time = iq_time[0:iq_idx_start + 1]
                                iq_dither = iq_dither[0:iq_idx_start + 1]
                                break
                            else:
                                excessive_time_count += 1
                        else:
                            excessive_time_count = 0
                    except:
                        pass #don't let one bad read bomb out

        else: #just the two nearest
            #now which to use? gc1 or gc2
            # base_tarfn = os.path.join(path, "gc1/gc1.tar")
            # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
            #     # try based on the virus_tar_path, up two levels
            #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc1.tar")

            base_tarfn = get_gc_path(cfg, "gc1", path)
            gc1_active=False
            if os.path.exists(base_tarfn):
                try:
                    #just checking for which is active
                    t1,p1 = Utils.open_file_from_tar(base_tarfn,gc1_near[0][0])
                    fh = fits.open(t1)
                    gc1_active = fh[0].header['GUIDLOOP'] == 'ACTIVE'
                    fh.close()
                    t1.close()
                except:
                    print(f"Exception in get_guider_fwhm: {traceback.format_exc()}")


            # base_tarfn = os.path.join(path, "gc2/gc2.tar")
            # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
            #     # try based on the virus_tar_path, up two levels
            #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc2.tar")

            base_tarfn = get_gc_path(cfg, "gc2", path)
            gc2_active=False
            if os.path.exists(base_tarfn):
                try:
                    #just checking for which is active
                    t1,p1 = Utils.open_file_from_tar(base_tarfn,gc2_near[0][0])
                    fh = fits.open(t1)
                    gc2_active = fh[0].header['GUIDLOOP'] == 'ACTIVE'
                    fh.close()
                    t1.close()
                except:
                    print(f"Exception in get_guider_fwhm: {traceback.format_exc()}")

            if gc1_active:
                # base_tarfn = os.path.join(path, "gc1/gc1.tar")
                # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
                #     # try based on the virus_tar_path, up two levels
                #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc1.tar")

                base_tarfn = get_gc_path(cfg, "gc1", path)
                gc_near = gc1_near

            elif gc2_active:
                # base_tarfn = os.path.join(path, "gc2/gc2.tar")
                # if not os.path.exists(base_tarfn) and cfg.virus_tar_path is not None:
                #     # try based on the virus_tar_path, up two levels
                #     base_tarfn = os.path.join(cfg.virus_tar_path.rstrip(".tar"), "gc2.tar")
                base_tarfn = get_gc_path(cfg, "gc2", path)
                gc_near = gc2_near
            else:
                print(f"No active guider. Cannot get seeing fwhm.")
                return None

            for near_list in gc_near:
                for near_file in near_list:
                    try:
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, near_file)
                        fh = fits.open(t1)
                        try:
                            iq.append(float(fh[0].header['IQ']))  # this is clipped later to exlcude 0s and other bad values
                        except:
                            try:
                                # card not there, try using PIXSCALE and OBFWHM
                                pixscale = float(fh[0].header['PIXSCALE'])
                                obfwhm = float(fh[0].header['OBFWHM'])
                                iq.append(obfwhm * pixscale)
                            except:
                                pass

                        fh.close()
                        t1.close()
                    except:
                        print(f"Invalid IQ card")

        if len(iq) == 1:
            np.savetxt(saved_fn,[iq[0]],fmt="%0.4f")
            return iq[0]
        elif len(iq) > 1:
            md_iq = np.nanmedian(np.clip(iq, 0.1, 9.0))
            np.savetxt(saved_fn, [md_iq], fmt="%0.4f")

            std_iq = np.nanstd(np.clip(iq, 0.1, 9.0))
            if len(iq_dither) > 0:
                with open(saved_detailed_fn,"w") as f:
                    f.write(f"Shot median: {md_iq:0.4f} +/- {std_iq:0.5f}\n")
                    for i in range(cfg.numexp):
                        sel = np.array(iq_dither) == i + 1
                        md_iq = np.nanmedian(np.clip(np.array(iq)[sel], 0.1, 9.0))
                        std_iq = np.nanstd(np.clip(np.array(iq)[sel], 0.1, 9.0))
                        f.write(f"  Dither #{i+1}: {md_iq:0.4f} +/- {std_iq:0.5f}\n")

                    f.write("All guider values:\n")
                    for i in range(len(iq_dither)):
                        f.write(f"{iq_dither[i]} {iq_time[i]} {iq[i]:0.4f}\n")

            return md_iq
        else:
            return None

    except:
        print(f"Exception in get_guider_fwhm: {traceback.format_exc()}")
        return None

def make_3d_friend_table_for_shot(detect_table, dsky_3D=6.0, dwave=4.0):
    """
    mostly lifted from hetdex_api

    using dwave to determine line vs cont (dwave == 0 ==> continuum source)

    :param detect_table:
    :param dsky_3D:
    :param dwave: (dwave == 0 ==> continuum source), note: all continuum sources get wave = 45
    :return:
    """


    if "line_detectid" in detect_table.columns:
        detid_col = "line_detectid"
        flux_col = "lineflux"
    elif "cont_detectid" in detect_table.columns:
        #this is assumed to be a continuum table or will ignore the wavelength
        if "wave" not in detect_table.columns:
            detect_table["wave"] = 4505.0 #set in the middle
        if "contflux" not in detect_table.columns:
            detect_table["contflux"] = [np.nanmedian(x[200:800]) for x in detect_table['obs_fluxd']]

        flux_col = "contflux"
        detid_col = "cont_detectid"
    elif "detectid" in detect_table.columns:
        detid_col = "detectid"
        if "contflux" in detect_table.columns:
            flux_col = "contflux"
        else:
            flux_col = "flux"
    else:
        print("ERROR! Unknown input table for make_3d_friend_table_for_shot()")
        return None


    kdtree, r = fof.mktree(
        detect_table["ra"],
        detect_table["dec"],
        detect_table["wave"],
        dsky=dsky_3D,
        dwave=dwave,
    )

    wfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(wfriend_lst) > 0:
        wfriend_table = fof.process_group_list(
            wfriend_lst,
            detect_table[detid_col],
            detect_table["ra"],
            detect_table["dec"],
            detect_table["wave"],
            detect_table[flux_col],
        )

        memberlist = []
        friendlist = []
        for row in wfriend_table:
            friendid = row["id"]
            members = np.array(row["members"])
            friendlist.extend(friendid * np.ones_like(members))
            memberlist.extend(members)

        wfriend_table.remove_column("members")

        wdetfriend_tab = Table()
        wdetfriend_tab.add_column(Column(np.array(friendlist), name="id"))
        wdetfriend_tab.add_column(Column(memberlist, name=detid_col))

        wdetfriend_shot = join(wdetfriend_tab, wfriend_table, keys="id")

        wdetfriend_shot.rename_column(detid_col, "detectid")
        wdetfriend_shot.rename_column("id", "wave_group_id")
        wdetfriend_shot.rename_column("size", "wave_group_size")
        wdetfriend_shot.rename_column("a", "wave_group_a")
        wdetfriend_shot.rename_column("b", "wave_group_b")
        wdetfriend_shot.rename_column("pa", "wave_group_pa")
        wdetfriend_shot.rename_column("icx", "wave_group_ra")
        wdetfriend_shot.rename_column("icy", "wave_group_dec")
        wdetfriend_shot.rename_column("icz", "wave_group_wave")

        wdetfriend_shot = wdetfriend_shot[
            "detectid",
            "wave_group_id",
            "wave_group_a",
            "wave_group_b",
            "wave_group_pa",
            "wave_group_ra",
            "wave_group_dec",
            "wave_group_wave",
        ]

        return wdetfriend_shot
    else:
        return None




def make_2d_friend_table_for_shot(detect_table,dsky_2D=3.0):
    """
    again, mostly lifted from hetdex_api
    :param detect_table:
    :return:
    """

    if "line_detectid" in detect_table.columns:
        detid_col = "line_detectid"
        flux_col = "lineflux"
    elif "cont_detectid" in detect_table.columns:
        #this is assumed to be a continuum table or will ignore the wavelength
        if "wave" not in detect_table.columns:
            detect_table["wave"] = 4505.0 #set in the middle
        if "contflux" not in detect_table.columns:
            detect_table["contflux"] = [np.nanmedian(x[200:800]) for x in detect_table['obs_fluxd']]

        flux_col = "contflux"
        detid_col = "cont_detectid"
    elif "detectid" in detect_table.columns:
        detid_col = "detectid"
        if "contflux" in detect_table.columns:
            flux_col = "contflux"
        else:
            flux_col = "flux"
    else:
        print("ERROR! Unknown input table for make_2d_friend_table_for_shot()")
        return None

    kdtree, r = fof.mktree(
        detect_table["ra"],
        detect_table["dec"],
        np.zeros_like(detect_table["ra"]),
        dsky=dsky_2D,
    )
    friend_lst = fof.frinds_of_friends(kdtree, r, Nmin=1)

    friend_table = fof.process_group_list(
        friend_lst,
        detect_table[detid_col],
        detect_table["ra"],
        detect_table["dec"],
        0.0 * detect_table["wave"],
        detect_table[flux_col],
    )

    memberlist = []
    friendlist = []
    for row in friend_table:
        friendid = row["id"]
        members = np.array(row["members"])
        friendlist.extend(friendid * np.ones_like(members))
        memberlist.extend(members)

    friend_table.remove_column("members")

    detfriend_tab = Table()
    detfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    detfriend_tab.add_column(Column(memberlist, name="detectid"))

    detfriend_shot = join(detfriend_tab, friend_table, keys="id")

    return detfriend_shot



def merge_wave_groups(tab, wid):

    try:
        sel_wid = tab["wave_group_id"] == wid
        grp = tab[sel_wid]
        sid, ns = np.unique(grp["source_id"], return_counts=True)

        sid_main = np.min(sid)

        sid_ind = []
        for sid_i in sid:
            for ind in np.where(tab["source_id"] == sid_i)[0]:
                sid_ind.append(ind)

        # now find any other wave groups and their associated source_id info to merge
        other_wids = np.unique(tab["wave_group_id"][sid_ind])

        for wid_i in other_wids:
            if wid_i == 0:
                continue
            elif wid_i == wid:
                continue

            sel_wid_i = tab["wave_group_id"] == wid_i
            grp = tab[sel_wid_i]
            sid, ns = np.unique(grp["source_id"], return_counts=True)

            for sid_i in sid:
                if sid_i == sid_main:
                    continue
                for ind in np.where(tab["source_id"] == sid_i)[0]:
                    sid_ind.append(ind)

        if np.size(np.unique(tab["source_id"][sid_ind])) > 1:
            return sid_ind
        else:
            return None
    except Exception:
        print("Merge wave group failed for {}".format(wid))
        return None


def precheck(cfg):
    """
    sanity check files, directories to see if this reduction can run
    (e.g. if the lib_calib path is not available for the date, then this cannot run)

    :param cfg:
    :return:
    """

    try:

        #if the user is (re)-running (e.g. --resume) and we are under a path that seems to be for
        # this shot, stop and warn. The thinking is that the user may have been checking on the reduction,
        # maybe editing progress.dat, and attempting to resume, but in the wrong directory
        try:
            if f"sci{cfg.datevshot}" in cfg.cwd_orig:
                print(f"[{cfg.datevshot}] WARNING! Will abort! Working directory may not be as intended: {cfg.cwd_orig}")
                print(f"[{cfg.datevshot}] It appears to already in include the reduction output dir (sci{cfg.datevshot}).")
                if cfg.resume:
                    print(f"[{cfg.datevshot}] Did you --resume from the wrong directory?")
                return -1
        except:
            print(f"Warning! Could not validate current working directory.",traceback.format_exc())


        # echo a few key paths:
        print(f"[{cfg.datevshot}] Precheck. HETDEX_API path: {hetdex_api_path}",flush=True)
        print(f"[{cfg.datevshot}] Precheck. ELiXer path: {elixer_path}",flush=True)

        month = cfg.datevshot[0:6]
        path_check = os.path.join(hetdex_projects_path,f"lib_calib/{month}")
        if not os.path.exists(path_check):
            print(f"[{cfg.datevshot}] Precheck fail. Dir does not exist: {path_check}")
            return -1

        if check_lib_calib(cfg) != 0:
            print(f"[{cfg.datevshot}] Precheck fail. lib_calib/{cfg.datevshot[0:6]} directory not ready")
            return -1

        #common missing installs (that don't show up until later)
        pkgs=['sklearn','numba']
        rc = 0
        for pkg in pkgs:
            if importlib.util.find_spec(pkg) is None:
                print(f"Fatal. You need to (pip) install '{pkg}'")
                rc = -1


        #warn if missing HTEDEX packages
        pkgs=["hetdex_shuffle","vdrp"]
        for pkg in pkgs:
            if importlib.util.find_spec(pkg) is None:
                print(f"Warn! You may need to manually install '{pkg}'")
                if pkg == "hetdex_shuffle":
                    install_help = """
                    for hetdex_shuffle (you need this BEFORE vdrp)
                        $ pip install --user svn+svn://luna.mpe.mpg.de/hetdexshuffle/trunk#egg=shuffle
                        then cd in that directory and run
                        $ pip install --user -e .
                    """
                    print(install_help)
                elif pkg == "vdrp":
                    install_help = """
                    for vdrp:
                        $ git clone https://github.com/HETDEX/vdrp.git
                        then cd in that directory and run
                        $ pip install --user -e .
                    """
                    print(install_help)

        #for hetdex_shuffle (you need this BEFORE vdrp)
        #$ pip install --user svn+svn://luna.mpe.mpg.de/hetdexshuffle/trunk#egg=shuffle
        #then cd in that directory and
        #$ pip install --user -e . (edited)


        #for vdrp:
        #$ git clone https://github.com/HETDEX/vdrp.git
        #then cd into it and run the
        #$ pip install --user -e .


        #is this a hetdex shot and if so, is it allowed?

        try:
            h5 = tables.open_file(HETDEXSurvey,mode="r")
            dex_shots = h5.root.Survey.read(field="shotid")
            h5.close()
            if np.int64(cfg.datevshot.replace("v","")) in dex_shots:
                cfg.hetdex_original = True
                if not cfg.hetdex and not cfg.multifits_only:
                    print(f"[{cfg.datevshot}] is an existing HETDEX shot. To re-reduce here, re-run with --hetdex")
                    rc = -1
                else:
                    #this is okay
                    print(f"[{cfg.datevshot}] is an existing HETDEX shot, but --hetdex or --multifits_only specified, so will re-reduce.")
        except:
            print(f"[{cfg.datevshot}] Exception checking HETDEX Survye file {HETDEXSurvey}", traceback.format_exc())


        if rc != 0:
            return -1

    except:
        print(f"Exception in precheck: {traceback.format_exc()}")

    return 0


def update_only(cfg):
    """

    :param cfg:
    :return:
    """

    lock = FileLock(Lock_mutex_fn)  # we are in the top directory (not sciXXXX)
    with lock:
        if cfg.update_local_repo:
            print("(Only) Updating local repo ... (this may take 1-2 minutes).")
            shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                            os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
            cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

            # if cfg.code_fn_to_copy is not None:
            #     if cfg.cwd_orig is not None:
            #         shutil.copy2(cfg.code_fn_to_copy,cfg.cwd_orig)
            #     else:
            #         shutil.copy2(cfg.code_fn_to_copy, os.getcwd())



def set_fitradecsp(cfg):
    try:
        os.environ['SSR_fitradecsp'] = f"{cfg.scriptdir}/fitradecsp"
        print(f"[{cfg.datevshot}] Set SSR_fitradecsp to {os.getenv('SSR_fitradecsp')}")
    except:
        print(f"Exception in set_fitradecsp: {traceback.format_exc()}")


def set_vred(cfg):
    try:
        os.environ['SSR_vred'] = f"{cfg.scriptdir}/vred"
        print(f"[{cfg.datevshot}] Set SSR_vred to {os.getenv('SSR_vred')}")
    except:
        print(f"Exception in SSR_vred: {traceback.format_exc()}")

def set_vred3(cfg):
    try:
        os.environ['SSR_vred3'] = f"{cfg.scriptdir}/vred3"
        print(f"[{cfg.datevshot}] Set SSR_vred3 to {os.getenv('SSR_vred3')}")
    except:
        print(f"Exception in SSR_vred3: {traceback.format_exc()}")

def get_het_raw_archive(cfg):
    """
    return best available het_raw archive
    prioritize /scratch over /work over /corral ... order should be set correctly
    in the global HETRaw_archive
    :return: path to the tar file (including the tarfile itself),
             True/False (True if is the datevshot.tar, False if date.tar)
             True/False (True if local (already copied), False if needs to be copied)
    """

    try:
        datevshot = cfg.datevshot

        #first, do we have it locally already? if so, it is always the datevshot.tar, not the higher level date.tar
        p = None
        if cfg.local_het_raw_path is None or len(cfg.local_het_raw_path)==0:
            cfg.local_het_raw_path = os.path.join(cfg.cwd_orig, "het_raw")

        if not os.path.exists(cfg.local_het_raw_path):
            Path(cfg.local_het_raw_path).mkdir(parents=True, exist_ok=True)

        p = os.path.join(cfg.local_het_raw_path,f"{datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")

        if p is not None:
            if os.path.exists(p):
                return p, True, True

        #otherwise, check the archival locations
        for p in HETRaw_archives:
            #check that it exists ... if so, stop
            if "/scratch" == p[:8]:
                #check for the path first
                fullpath = os.path.join(p,f"{datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                if os.path.exists(fullpath):
                    return fullpath,True,False

                #might just be the day tar file
                fullpath = os.path.join(p,f"{datevshot[:8]}.tar")
                if os.path.exists(fullpath):
                    return fullpath,False,False

            elif "/work" == p[:5]: #this is the usual, but incomplete spot and might have a tar or an exploded path
                #check for the path first
                fullpath = os.path.join(p,f"{datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                if os.path.exists(fullpath):
                    return fullpath,True,False

                #might just be the day tar file
                fullpath = os.path.join(p,f"{datevshot[:8]}.tar")
                if os.path.exists(fullpath):
                    return fullpath,False,False
            elif "/corral" == p[:7]:  #notice includes /corral and /corral-repl
                #only day tars
                fullpath = os.path.join(p,f"{datevshot[:8]}.tar")
                if os.path.exists(fullpath):
                    return fullpath,False,False
            else: #this is a problem
                print(f"[{cfg.datevshot}] Problem identifying het_raw archive. Not found.")
                return None, False, False

        print(f"[{cfg.datevshot}] Problem identifying het_raw archive. Could not locate archive.")
        return None, False, False #did not find the file anywhere ... this is an error
    except:
        print(f"[{cfg.datevshot}] Exception in get_het_raw_archive: {traceback.format_exc()}")
        return None, False, False

def get_het_raw_copy_mutex(cfg,virus_path,release=None):
    """

                # we need to limit our access hits to corral across ALL instances
                # so there is an additional mutex
                # the following waits on a mutex to edit a counter
                # if the counter is okay, update to include THIS copy
                #    and release the mutex
                #if not okay, releast the mutex, sleep then try again
                #Once copy is complete, wait on mutex, decriment the counter

    :param cfg:
    :param virus_path:
    :param release: decrement this counter file or if, None, obtain a lock and increment if needed
    :return: the counter file to decrement (or None) if not needed
    """

    global CorralMaxCopyLimit, WorkMaxCopyLimit

    try:
        if "/corral" == virus_path[:7]: #needs a limit
            mux_root = "/".join(cfg.cwd_orig.split("/")[:4]) + "/ssr_mux/corral"
            if not os.path.exists(mux_root):
                Path(mux_root).mkdir(parents=True, exist_ok=True)

            if release is None: #starting a copy request
                #need a user level storage location for the mutex and the counter so it does not matter
                #how many nodes or jobs are in use, they will all check here
                # I think I can assume the top level to be the first three tokens of the cwd()

                io_lock_file = os.path.join(mux_root,"ssr_io.lock")
                redlight = True

                # no PID here and no shot, just date ... want it common to any of the save datevshot
                sync_fn = f"{cfg.datevshot[:8]}.sync"
                while redlight:
                    lock = FileLock(io_lock_file)
                    with lock:
                        fns = glob.glob(os.path.join(mux_root, "*.sync"))
                        active = len(fns)
                        if active < CorralMaxCopyLimit: #strictly less than since THIS one will be +1
                            #what if this file is already being copied/untarred ... also need to wait
                            #since this same file is date based, multiple datevshots could be trying to hit it
                            #at the same time
                            if os.path.exists(os.path.join(mux_root, sync_fn)):
                                line = None
                                try:
                                    with open(os.path.join(mux_root, sync_fn), "r") as f:
                                        line = f.readline()
                                except:
                                    line = "Unable to read .sync file"
                                print(f"[{cfg.datevshot}] copy must wait. Archive tar already being copied by another process: {line}")
                            else:
                                with open(os.path.join(mux_root,sync_fn), "w") as f:
                                    f.write(f"PID: {os.getpid()} BEGAN {str(datetime.now())} from {cfg.cwd_orig}\n")
                                    print(f"[{cfg.datevshot}] cleared to copy from /corral. {active} other active copies.")
                                    redlight = False
                        else:
                            print(f"[{cfg.datevshot}] too many active ({active}) copies from /corral. Limit ({CorralMaxCopyLimit}) . Must wait ...")

                    if redlight:
                        time.sleep(SafeActiveShotsSleep)

                return sync_fn #name of the sync file
            else:  # we are done with the copy, just delete our sync file
                try:
                    os.remove(os.path.join(mux_root, release))
                except:
                    print(f"[{cfg.datevshot}] Exception in get_het_raw_copy_mutex, release corral: {traceback.format_exc()}")
        elif "/work" == virus_path[:5]:  # needs a limit
            mux_root = "/".join(cfg.cwd_orig.split("/")[:4]) + "/ssr_mux/work"
            if not os.path.exists(mux_root):
                Path(mux_root).mkdir(parents=True, exist_ok=True)

            if release is None:  # starting a copy request
                # need a user level storage location for the mutex and the counter so it does not matter
                # how many nodes or jobs are in use, they will all check here
                # I think I can assume the top level to be the first three tokens of the cwd()

                io_lock_file = os.path.join(mux_root, "ssr_io.lock")
                redlight = True

                #these are datevshot based, not just date
                sync_fn = f"{cfg.datevshot}.sync" #no PID here, want it common to any of the same datevshot
                while redlight:
                    lock = FileLock(io_lock_file)
                    with lock:
                        fns = glob.glob(os.path.join(mux_root, "*.sync"))
                        active = len(fns)
                        if active < WorkMaxCopyLimit:  # strictly less than since THIS one will be +1

                            #what if this file is already being copied ... also need to wait
                            if os.path.exists(os.path.join(mux_root, sync_fn)):
                                line = None
                                try:
                                    with open(os.path.join(mux_root, sync_fn), "r") as f:
                                        line = f.readline()
                                except:
                                    line = "Unable to read .sync file"
                                print(f"[{cfg.datevshot}] copy must wait. Archive tar already being copied by another process: {line}")
                            else:
                                with open(os.path.join(mux_root, sync_fn), "w") as f:
                                    f.write(f"PID: {os.getpid()} BEGAN {str(datetime.now())} from {cfg.cwd_orig}\n")
                                    print(
                                        f"[{cfg.datevshot}] cleared to copy from /work. {active} other active copies.")
                                    redlight = False
                        else:
                            print(
                                f"[{cfg.datevshot}] too many active ({active}) copies from /work. Limit ({WorkMaxCopyLimit}) . Must wait ...")

                    if redlight:
                        time.sleep(SafeActiveShotsSleep)

                return sync_fn  # name of the sync file
            else: #we are done with the copy, just delete our sync file
                try:
                    os.remove(os.path.join(mux_root, release))
                    #os.remove(release)
                except:
                    print(f"[{cfg.datevshot}] Exception in get_het_raw_copy_mutex, release work: {traceback.format_exc()}")
        else:
            return None #no limit needed
    except:
        print(f"[{cfg.datevshot}] Exception in get_het_raw_copy_mutex: {traceback.format_exc()}")



def copy_het_raw_file(cfg):
    """
    handle the copying of the het raw file with mutex and datevshot.tar vs date.tar
    also used to wait on a copy that is in progress
    :param cfg:
    :return:
    """

    rc = 0
    try:
        virus_path, is_dvs_tar, is_local = get_het_raw_archive(cfg)
        if virus_path is None: #problem, cannot continue
            print(f"[{cfg.datevshot}] Fatal! Cannot find het_raw data.")
            return -1


        if is_local: #already have it OR it is being copied ... check for the mutex

            cfg.virus_tar_path = virus_path
            #if it is being copied or untarred, the lock will be held by the other process that is doing it
            #so, wait on the lock ... once we have it, the copy or untar should be complete
            cfg.copy_lock_file = os.path.join(cfg.local_het_raw_path, f"{cfg.datevshot}.lock")
            lock = FileLock(cfg.copy_lock_file)
            with lock:  # we have the lock now
                destination_path = os.path.join(cfg.local_het_raw_path,
                                                f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                if os.path.exists(destination_path):  # good copy
                    print(f"[{cfg.datevshot}] Using {destination_path}")
                else:
                    print(f"[{cfg.datevshot}] Unexpected error. {destination_path} does not exist")
                    return -1

            # copy_check = virus_path[:-3] + "copy"  #always ends in .tar, use .copy to indicate in progress copy
            # #use this as both the mutex and the copy in progress indicator
            # while os.path.exists(copy_check):
            #     #some other process is copying it ... wait for it to be done.
            #     print(f"[{cfg.datevshot}] Waiting on copy of: {virus_path}")
            #     time.sleep(SafeActiveShotsSleep) #reuse the sleep delay
            # print(f"[{cfg.datevshot}] Using {virus_path}")
            return 0 #ready to go
        else: #this needs to be copied
            #get the copy mutex ... this also handles limits on /corral

            #once we have the mutex ... check again if the file already exists ... could be another
            #process got to it first and it is ready now


            #copy_lock_file = os.path.join(cfg.local_het_raw_path, os.path.basename(virus_path)[0:-3] +"lock")
            cfg.copy_lock_file = os.path.join(cfg.local_het_raw_path,  f"{cfg.datevshot}.lock")
            lock = FileLock(cfg.copy_lock_file)
            with lock: #we have the lock now
                #does the tar file already exist? someone else already copied got it?
                # if os.path.exists(virus_path):  # someone else already copied it.
                #     print(f"[{cfg.datevshot}] Using (already copied) {virus_path}")

                # we need to limit our access hits to corral across ALL instances
                # so there is an additional mutex
                # the following waits on a mutex to edit a counter
                # if the counter is okay, update to include THIS copy
                #    and release the mutex
                #if not okay, releast the mutex, sleep then try again
                #Once copy is complete, wait on mutex, decriment the counter

                destination_path = os.path.join(cfg.local_het_raw_path,
                                                f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                if os.path.exists(destination_path):  # good copy ... someone else already handled it
                    print(f"[{cfg.datevshot}] Using {destination_path}")
                    cfg.virus_tar_path = destination_path
                else: #do the copy or untar

                    counter_file = get_het_raw_copy_mutex(cfg, virus_path)

                    #start the copy
                    #make sure the recieving location exists
                    if not os.path.exists(cfg.local_het_raw_path):
                        os.makedirs(cfg.local_het_raw_path, exist_ok=True)

                    # if True:
                    #     ct = 0
                    #     while ct < 10:
                    #         ct += 1
                    #         destination_path = os.path.join(cfg.local_het_raw_path,
                    #                                         f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                    #         print(f"*** DEBUG FAKE COPY *** {virus_path} to {destination_path}")
                    #         time.sleep(30.0)
                    # else:
                    if is_dvs_tar:
                        # IF this is the datevshottar then just copy it (probably this is on /work)
                        # still need to check for the gc1 and gc2 stuff
                        destination_path = os.path.join(cfg.local_het_raw_path,
                                                        f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")
                        Path(os.path.dirname(destination_path)).mkdir(parents=True, exist_ok=True)
                        cmd = f"cp {virus_path} {destination_path}"
                        system_command(cfg, cmd)
                        if os.path.exists(destination_path): #good copy
                            cfg.virus_tar_path = destination_path
                            try:
                                ix = virus_path.index(cfg.datevshot[:8])
                                basepath = virus_path[:ix] + f"{cfg.datevshot[:8]}/"
                            except:
                                basepath = None

                            if basepath is not None:
                                destination_path = os.path.join(cfg.local_het_raw_path, f"{cfg.datevshot[:8]}/gc1")
                                Path(destination_path).mkdir(parents=True, exist_ok=True)

                                #!! since the gc?.tar (or images.tar) are part of the top level date
                                #   and not tied directly to datevshot, they might already exist locally even if this
                                #   specific datevshot did not
                                if len(glob.glob(destination_path + "/*.tar")) == 0:
                                    #e.g. /work/03946/hetdex/maverick/20260411/gc1
                                    # NOTICE: *.tar as sometimes it is gc2.tar or images.tar
                                    source_path = os.path.join(basepath,"gc1/*.tar")
                                    cmd = f"cp {source_path} {destination_path}"
                                    system_command(cfg, cmd)

                                destination_path = os.path.join(cfg.local_het_raw_path, f"{cfg.datevshot[:8]}/gc2")
                                Path(destination_path).mkdir(parents=True, exist_ok=True)
                                if len(glob.glob(destination_path + "/*.tar")) == 0:
                                    #e.g. /work/03946/hetdex/maverick/20260411/gc2
                                    #NOTICE: *.tar as sometimes it is gc2.tar or images.tar
                                    source_path = os.path.join(basepath,"gc2/*.tar")
                                    cmd = f"cp {source_path} {destination_path}"
                                    system_command(cfg, cmd)

                        else:
                            print(f"[{cfg.datevshot}] (1) Failed to copy/extract date tar file.", traceback.format_exc())
                            rc = -1
                    else:
                        #this is the date.tar, then need to extract parts of it
                        try:
                            destination_path = cfg.local_het_raw_path
                            file_to_extract = f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar"
                            Path(destination_path).mkdir(parents=True, exist_ok=True)
                            cmd = f"tar -xvf {virus_path} -C {destination_path} {file_to_extract}"
                            system_command(cfg, cmd)
                            if os.path.exists(destination_path): #good copy, continue with the other 2 but don't check them
                                cfg.virus_tar_path = destination_path
                                # this MIGHT already exist if another process got it
                                #  since the gc?.tar (or images.tar) are part of the top level date
                                #  and not tied directly to datevshot, they might already exist locally even if this
                                #  specific datevshot did not
                                #destination_path = os.path.join(cfg.local_het_raw_path,f"{cfg.datevshot[:8]}/gc1")
                                file_to_extract = f"{cfg.datevshot[:8]}/gc1/*.tar"
                                #Path(destination_path).mkdir(parents=True, exist_ok=True)
                                if not os.path.exists(destination_path) or len(glob.glob(destination_path + "/*.tar")) == 0:
                                    cmd = f"tar -xvf {virus_path} -C {destination_path} {file_to_extract}"
                                    system_command(cfg, cmd)

                                #destination_path = os.path.join(cfg.local_het_raw_path,f"{cfg.datevshot[:8]}/gc2")
                                file_to_extract = f"{cfg.datevshot[:8]}/gc2/*.tar"
                                #Path(destination_path).mkdir(parents=True, exist_ok=True)
                                if not os.path.exists(destination_path) or len(glob.glob(destination_path + "/*.tar")) == 0:
                                    cmd = f"tar -xvf {virus_path} -C {destination_path} {file_to_extract}"
                                    system_command(cfg, cmd)
                            else:
                                print(f"[{cfg.datevshot}] (2) Failed to copy/extract date tar file.", traceback.format_exc())
                                rc = -1
                        except:
                            print(f"[{cfg.datevshot}] (3) Failed to copy/extract date tar file.", traceback.format_exc())
                            rc = -1

                    #release the mutex/counter decrement
                    if counter_file is not None:
                        get_het_raw_copy_mutex(cfg, virus_path,release=counter_file)

                # if True:
                #     print(f"*** DEBUG FORCE EXIT *** ")
                #     rc = -1
    except:
        print(f"[{cfg.datevshot}] Exception in copy_het_raw_file: {traceback.format_exc()}")
        rc = -1

    return rc

def initial_setup(cfg):
    """
    copy from script repo(s)

    change the cwd to the new workdir for this shot

    this is largely equivalent to what rsetups would do, but here for a single shot and not a whole month

    :return:
    """

    workdir = os.path.join(WorkDirRoot,f"sci{cfg.datevshot}")

    repair = False
    resume = False #notice: cfg.resume MAY be true and this can still be false if the directory does not already exist
    if os.path.exists(workdir):
        if cfg.repair:
            print(f"[{cfg.datevshot}] Resume + Repair. Leave directory intact, but re-copy work files: {workdir}")
            repair = True
        elif cfg.resume:
            print(f"[{cfg.datevshot}] Resuming. Leave directory intact: {workdir}")
            resume = True
        elif cfg.overwrite:
            print(f"[{cfg.datevshot}] Overwriting directory {workdir} ... ")
            shutil.rmtree(workdir)
        else:
            print(f"[{cfg.datevshot}] Shot directory already exists here! {workdir}")
            print(f"[{cfg.datevshot}] Please include --resume or --overwrite to make intention clear.")
            print(f"[{cfg.datevshot}] Or is this a repeated SLURM task?")
            if cfg.clean > 1:
                print(f"[{cfg.datevshot}] Did you intend to only clean up the directory? "
                      f"If so, --clean needs to be the negative value. (see --help)")
            return -1

    if not resume:
        #os.makedirs(workdir)
        Path(workdir).mkdir(parents=True, exist_ok=True)

        if LocalScriptRepo is not None:
            #need to see if can get the file lock in case another instance is copying
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:
                if cfg.update_local_repo:
                    print("Updating local repo ... (this may take 1-2 minutes)")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                    # if cfg.code_fn_to_copy is not None:
                    #     if cfg.cwd_orig is not None:
                    #         shutil.copy2(cfg.code_fn_to_copy, cfg.cwd_orig)
                    #     else:
                    #         shutil.copy2(cfg.code_fn_to_copy, os.getcwd())

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print(f"[{cfg.datevshot}] Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    #copy first to local script repo
                    print(f"[{cfg.datevshot}] Copying to local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(),LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

            #lock auto releases
        else:
            print(f"[{cfg.datevshot}] Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")
    elif repair:
        if LocalScriptRepo is not None:
            fatal_rtn = False
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:
                if cfg.update_local_repo:
                    print(f"[{cfg.datevshot}] Updating local repo ... (this may take 1-2 minutes)")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print(f"[{cfg.datevshot}] Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    print(f"[{cfg.datevshot}] Fatal! --repair selected, but no script repo.")
                    fatal_rtn = True
            # lock auto releases
            if fatal_rtn:
                return -1
        else:
            print(f"[{cfg.datevshot}] Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")
    else:
        if LocalScriptRepo is not None:
            fatal_rtn = False
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:

                if cfg.update_local_repo:
                    print(f"[{cfg.datevshot}] Updating local repo ... (this may take 1-2 minutes)")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print(f"[{cfg.datevshot}] Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    print(f"[{cfg.datevshot}] Fatal! --resume selected, but no script repo.")
                    fatal_rtn = True
            # lock auto releases
            if fatal_rtn:
                return -1
        else:
            print(f"[{cfg.datevshot}] Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")

    set_fitradecsp(cfg)
    set_vred(cfg)
    set_vred3(cfg)
    os.chdir(workdir)
    cfg.cwd = os.getcwd() #now under the sci<shot> directory

    #check for virus.tar files

    #if True: #new way
    if copy_het_raw_file(cfg) < 0:
        #FATAL
        print(f"[{cfg.datevshot}] FATAL! Could not obtain het_raw data.")
        return -1
    # else: #old way
    #     virustar = f"{cfg.datevshot[0:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar"
    #     if cfg.local_het_raw_path is None or len(cfg.local_het_raw_path) == 0:
    #         cfg.local_het_raw_path = os.path.join(cfg.cwd_orig,"het_raw")
    #     else:
    #         if os.path.basename(cfg.local_het_raw_path) != "het_raw":
    #             cfg.local_het_raw_path = os.path.join(cfg.local_het_raw_path,"het_raw")
    #     virus_paths = [HET_by_date,cfg.local_het_raw_path] #, os.path.join(cfg.cwd_orig,"het_raw")]
    #     if not os.path.exists(os.path.join(virus_paths[0],virustar)):
    #         if not os.path.exists(os.path.join(virus_paths[1],virustar)):
    #             fail = True
    #             src_tar = os.path.join(HETRaw_archive,f"{cfg.datevshot[:8]}.tar")
    #             if os.path.exists(src_tar):
    #                 print(f"[{cfg.datevshot}]. Could not locate {virustar} under {virus_paths}")
    #                 print(f"[{cfg.datevshot}]. Will attempt to copy {src_tar}. If successful, this will take many minutes ...")
    #
    #                 #tmp lock file (so only one attempt to copy and extract; since there are often many observations
    #                 #   in this file, only one of the tasks that are on that date should copy)
    #                 tmp_lock_file =f"{cfg.datevshot[:-4]}.lock"
    #                 lock = FileLock(f"{cfg.cwd_orig}/{tmp_lock_file}")  # we are in sciXXX, but want the lock file up a level
    #                 with lock:
    #                     #try to copy
    #                     #if not os.path.exists(f"{cfg.cwd_orig}/het_raw"):
    #                     #    os.makedirs(f"{cfg.cwd_orig}/het_raw", exist_ok=True)
    #                     if not os.path.exists(cfg.local_het_raw_path):
    #                         os.makedirs(cfg.local_het_raw_path, exist_ok=True)
    #
    #                     if safe_cd(cfg.local_het_raw_path):
    #                         #check that another process has not already copied the file (see above pathing)
    #                         if not os.path.exists(os.path.join(virus_paths[1], virustar)):
    #
    #                             #shutil.copy2(os.path.join(HETRaw_archive,f"{cfg.datevshot[:8]}.tar"),".")
    #
    #                             try:
    #                                 # need only some of the paths ... don't know which /virus we may also need so keep all
    #                                 #cmd = f"tar -xvf {cfg.datevshot[:8]}.tar {cfg.datevshot[:8]}/virus"
    #                                 cmd = f"tar -xvf {src_tar} {cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar"
    #                                 system_command(cfg, cmd)
    #
    #                                 #cmd = f"tar -xvf {cfg.datevshot[:8]}.tar gc1"
    #                                 #this MIGHT already exist if another process got it
    #                                 if not os.path.exists(f"{cfg.datevshot[:8]}/gc1"):
    #                                     cmd = f"tar -xvf {src_tar} {cfg.datevshot[:8]}/gc1"
    #                                     system_command(cfg, cmd)
    #
    #                                 #cmd = f"tar -xvf {cfg.datevshot[:8]}.tar gc2"
    #                                 if not os.path.exists(f"{cfg.datevshot[:8]}/gc2"):
    #                                     cmd = f"tar -xvf {src_tar} {cfg.datevshot[:8]}/gc2"
    #                                     system_command(cfg, cmd)
    #
    #                                 #last check (note: gc1 and gc2 are desirable to have, but not required)
    #                                 if os.path.exists(os.path.join(virus_paths[1], virustar)):
    #                                     fail = False #we should be okay now
    #
    #                             except:
    #                                 print(f"[{cfg.datevshot}] Failed to copy/extract date tar file.",traceback.format_exc())
    #
    #                             #not necessary ... changed to just extract the key files
    #                             #now we should delete the big <date>.tar file
    #                             #cmd = f"rm  {cfg.datevshot[:8]}.tar"
    #                             #system_command(cfg, cmd)
    #
    #                         else:
    #                             fail = False #another process DID copy and the file we want is there now
    #
    #             if fail:
    #                 print(f"[{cfg.datevshot}] FATAL. Could not locate {virustar} under {virus_paths}")
    #                 print(f"[{cfg.datevshot}] You may need to first copy and extract "
    #                       f"{HETRaw_archive}/<date>.tar to your local het_raw directory")
    #                 return -1
    #             else: #we did eventually get what we need, so go back to the working dir and continue
    #                 os.chdir(cfg.cwd)
    #                 cfg.virus_tar_path = os.path.join(cfg.local_het_raw_path,virustar)
    #         else:
    #             cfg.virus_tar_path = os.path.join(virus_paths[1], virustar)
    #     else:
    #         cfg.virus_tar_path = os.path.join(virus_paths[0],virustar)

    if not resume or cfg.update_local_repo or repair:

        # some other process might be updating the local repo ... do not proceed IF there is a lock
        lock = FileLock(Lock_mutex_fn)  # we are in the top directory (not sciXXXX)
        with lock:
            # obtained the lock, so we should be good now,
            # just release and go
            print(f"[{cfg.datevshot}] Mutex checked. Okay. Safe to local_repo.")
            # lock auto releases

        print(f"[{cfg.datevshot}] Copying source code to working directory {cfg.cwd}...")
        ## if ANY of this fails it is fatal

        #shutil.copy2(os.path.join(cfg.scriptdir, "science_reductions", "rsetups"),".") #no, this function is its equivalent

        shutil.copy2(os.path.join(cfg.scriptdir, "rfixspec"), ".")
        shutil.copytree(os.path.join(cfg.scriptdir, "sciscripts"), ".", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir,"vdrp"), "vdrp", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "detect"), "detect", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "getcen"), "getcen", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "alldet"), "alldet", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "cs"), "cs", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "Diagnose"), "Diagnose", dirs_exist_ok=True)

        #update the "home" path tilde
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits")
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_fix")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_field")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_fix")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbacks")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbacks_2")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_s")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rimarb")  # use '#' as sed separator rather than "/"


        #update the old red1 path; should now all point to the common directory above each sci<YYYYMMDDvSSS> directories
        #yes, I want cwd_orig ... I want a single "reductions" directory off the top with all the sciXXX as siblings
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_field")  # use '#' as sed separator rather than "/"
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_fix")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rerun2")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rtaremc") #not necessary, just for completeness
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtar")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtarm.defunct") #not necessary, just for completeness
        #next two are just safeties, the scripts THIS code is using should already have them updated
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# alldet/rfft")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# getcen/rgetifucen")

        #extra files needed
        #vdrp : /work/00115/gebhardt/maverick/fplane ... need the fplane for the date
        os.makedirs(os.path.join(cfg.cwd,"vdrp/fplane"), exist_ok=True)


        if os.path.exists(os.path.join(karlfplane, f"fp{cfg.datevshot[0:8]}")):
            shutil.copy2(os.path.join(karlfplane, f"fp{cfg.datevshot[0:8]}"), os.path.join(cfg.cwd,"vdrp/fplane"))
        else:
            print(f"[{cfg.datevshot}] !Warning! fplane file (fp{cfg.datevshot[0:8]}) not found. Using last known ({LastKnownFplane}) instead. ")
            shutil.copy2(os.path.join(karlfplane, f"{LastKnownFplane}"), os.path.join(cfg.cwd, f"vdrp/fplane/fp{cfg.datevshot[0:8]}"))

        #fix paths in the . cfg files
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.original")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.gaia")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.sdss")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.panstarrs")

        system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config")
        system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.original")
        system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.gaia")
        system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.sdss")
        system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.panstarrs")

        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.original")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.gaia")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.sdss")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.panstarrs")

        system_command(cfg, f"sed -i s#/scratch/00115/gebhardt#{cfg.cwd}# vdrp/shifts/runsh1")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/shifts/runsh2")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# vdrp/shifts/run_shifts.sh")

        #vdrp: need the runshifts for the shot (from karlgettar)
        #e.g. run_shifts.sh 20240730 009 16.317927 33.689304 1
        # (formerly, this was the "clean_rta" script

        if os.path.exists(os.path.join(karlgettar,f"rta.{cfg.datevshot[0:6]}")):
            #rta_v is the track 0 or 1 for east/ west?
           rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(os.path.join(karlgettar,f"rta.{cfg.datevshot[0:6]}"),
                                                                usecols=[1,2,3,4,5],unpack=True,dtype=str)
        else:
            date_year = int(cfg.datevshot[0:4])
            #if 202407 < int(cfg.datevshot[0:6]) <= 202412:
            if date_year == 2024:
                rtafn = os.path.join(karlgettar,f"rta.202488")
            elif date_year >= 2025:
                rtafn = os.path.join(karlgettar, f"rta.{date_year}00")
                #rtafn = os.path.join(karlgettar, f"rta.202500")
            else:
                #should not happen
                rtafn = None #just to turn off warning
                Quit(cfg,-1,"Fatal. Cannot locate suitable rta file")

            rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(rtafn,
                                                                usecols=[1, 2, 3, 4, 5], unpack=True, dtype=str)

        sel = (rta_date == cfg.datevshot[0:8]) * (rta_shot == cfg.datevshot[-3:])
        if np.count_nonzero(sel) == 0:
            #none found, so make up my own .... but we don't have this info yet
            cfg.build_rta = True


        else:
            with open(os.path.join(cfg.cwd,f"vdrp/shifts/rta.{cfg.datevshot[0:6]}"),"w") as f:
                for d,s,ra,dec,v in zip(rta_date[sel], rta_shot[sel], rta_ra[sel], rta_dec[sel], rta_v[sel]):
                    #note ra here is in decimal hours
                    f.write(f"run_shifts.sh {d} {s} {ra} {dec} {v} \n")


        #fix detect stuff
        system_command(cfg, f"sed -i s/rsetstar/#/ detect/rallcal") #comment out rsetstar call as we want to run it separately
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rfitfw0")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsetstar")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsp3f")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsp3fc")


        if cfg.repair: #from here on out, resume and repair are the same
            cfg.resume = True
        return 0
    #end if not resume

    return 0



def add_text_to_image(cfg, image_path: str, text: str):
    """
    Opens a PNG image, adds text to the upper right corner, and saves it.

    Args:
        image_path: Path to the input PNG file
        text: Text to add to the image
        output_path: Path to save the output (defaults to overwriting the input)
    """


    try:
        # Open the image
        img = Image.open(image_path).convert("RGBA")
        draw = ImageDraw.Draw(img)

        # Try to use a truetype font, fall back to default if not available
        font_size = max(20, img.width // 30)  # Scale font size to image width
        try:
            #this is the default path on linux systems
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", font_size)
        except IOError:
            #this will not work, does not give you a TrueType font needed to insert into png
            #font = ImageFont.load_default()
            print(f"[{cfg.datevshot}] Exception! in add_text_to_image(). Cannot find fonts.", traceback.format_exc())
            return

        # Measure text size
        bbox = draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]

        # Position: upper right corner with a small margin
        margin = 10
        x = img.width - text_width - margin
        y = margin

        # Draw a semi-transparent background rectangle for readability
        padding = 5
        bg_bbox = (x - padding, y - padding, x + text_width + padding, y + text_height + padding)
        overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
        overlay_draw = ImageDraw.Draw(overlay)
        overlay_draw.rectangle(bg_bbox, fill=(0, 0, 0, 140))
        img = Image.alpha_composite(img, overlay)

        # Draw the text in white
        draw = ImageDraw.Draw(img)
        draw.text((x, y), text, font=font, fill=(255, 255, 255, 255))


        # Save — convert back to RGBA-safe PNG
        img.save(image_path, format="PNG")
        print(f"[{cfg.datevshot}] Updated image: {image_path}")
    except:
        print(f"[{cfg.datevshot}] Exception! in add_text_to_image()", traceback.format_exc())


def get_avg_sky(cfg):
    """
    can run AFTER run1s is done
    looks at sciXXXX/alldet/dXXXamp.dat
    and gets averge of all IFUs for columns Avg_orig

    values > 1000.0 are a problem
    values >= 9999.0 are fatal


    note: the d*amp.dat file under alldet/output is copied from /tmp (amp.out) during rfft run,
             in the overall step4 ... so, fairly late, AFTER vdrp and not long before the h5 is created
             *** you cannot stop at --multifits_only and get this
    :param cfg:
    :return:
    """
    avg_sky = None
    try:
        #not sure where I want to call this, so lets use the top level as a base path
        #filename like:  d20250101s024exp01amp.dat
        pattern = os.path.join(cfg.cwd_orig,f"sci{cfg.datevshot}/alldet/output/d*amp.dat")
        fns = glob.glob(pattern)

        #filename like: d20250101s024exp011de.png
        pattern = os.path.join(cfg.cwd_orig, f"sci{cfg.datevshot}/d{cfg.datevshot[0:8]}*/d*.png")
        img_fns = glob.glob(pattern)
        img_exps = [os.path.basename(x).split("exp")[1][:2] for x in img_fns]

        if len(fns) == 0:
            print(f"[{cfg.datevshot}] Can not collect avg_sky. Could not locate: {pattern}")
            return None
        #column headers
        #Spc_slt_iid_am Factor N_c Avg Scale W0 W1 Nlo Avg_orig chi Frac_c2 Frac0
        avg = []

        for fn in fns:
            x_avg = np.nanmean(np.loadtxt(fn,dtype=float,usecols=8,skiprows=1,unpack=True))
            if x_avg is not None:
                avg.append(x_avg)
                try:
                    #which exposure:
                    exp = os.path.basename(fn).split('exp')[1][:2] #should be a 2 digit string
                    ix = img_exps.index(exp)
                    img_path = img_fns[ix]
                    add_text_to_image(cfg, img_path, f"sky {x_avg:0.1f}")
                except:
                    print(f"[{cfg.datevshot}] Exception! in get_avg_sky(), trying to add sky value to png", traceback.format_exc())

        if len(avg) > 1:
            avg_sky = np.nanmean(avg)
        elif len(avg) == 1:
            avg_sky = avg[0]
        else: #this is a problem
            print(f"[{cfg.datevshot}] Warning! Could not compute an average sky in get_avg_sky(): avg={avg}")

        #the lower level file may be deleted in the clean step, so save off here so can
        #update the h5 file in the ssr_hdf5 call in the next process
        if avg_sky is not None:
            with open(os.path.join(cfg.cwd_orig,f"sci{cfg.datevshot}/avg_sky.dat"),"w") as f:
                f.write(f"{avg_sky}\n")

    except:
        print(f"[{cfg.datevshot}] Exception! in get_avg_sky()", traceback.format_exc())

    return avg_sky

def get_detection_counts(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        fn = os.path.join(cfg.cwd_orig,f"sci{cfg.datevshot}/{cfg.datevshot}_line_sourcecat.tab")
        if os.path.exists(fn):
            t = Table.read(fn,format="ascii")
            if t is not None:
                cfg.num_line_dets = len(t)
                del t
        else:
            print(f"[{cfg.datevshot}] Warning! Could not load number of line detections. Not found: {fn}")

        fn = os.path.join(cfg.cwd_orig, f"sci{cfg.datevshot}/{cfg.datevshot}_cont_sourcecat.tab")
        if os.path.exists(fn):
            t = Table.read(fn,format="ascii")
            if t is not None:
                cfg.num_cont_dets = len(t)
                del t
        else:
            print(f"[{cfg.datevshot}] Warning! Could not load number of cont detections. Not found: {fn}")
    except:
        print(f"[{cfg.datevshot}] Warning! Could not load number of detections.",traceback.format_exc())




def update_vdrp_config_limits(cfg):
    """
    update the vmin/vmax stretch on the collapsed imaging and the faint maglimit for calibration stars
    based on the total exposure time

    the concern is that, with long exposures we get too many faint stars and it can get confusing, so we cap the limit

    :param cfg:
    :return:
    """

    try:
        if cfg.total_exp_time is not None and cfg.total_exp_time > 1800.0:

            # if cfg.total_exp_time is not None:
            #     numexp = max(1, cfg.numexp)
            #     multiplier = max(1.0, (cfg.total_exp_time / numexp / 360.))  # 360secs is the nominal default for a HETDEX exposure
            #     tp_upper = max(0.25, 0.25 * multiplier)
            # else:
            #     multiplier = 1.0

            multiplier = linear_exptime_scale(cfg)

            #figure out what to use: should go roughly as sqr of time, but lets go linear here
            #multiplier = (cfg.total_exp_time / 1080.0) ** 2 #HETDEX 3Dither standard
            magmax = round(22.0 - 2.5 * np.log10(multiplier), 1) #22.0 is mktot_magmax (note 0.0 is always the min)
            #range is 50 so, centered at 5 (+/-25)
            # vmin = round(-20. * multiplier**2) #-20 is cofes_vis_vmin default
            # vmax = round(30. * 100) #multiplier**2) #30 is cofes_vis_vmax default
            # ctr = round(5.0 * multiplier)
            # vmin = ctr - 25 * multiplier
            # vmax = ctr + 25 * multiplier

            #vmin = 0 #if both are zero, it will trigger a zScale calculation
            #vmax = 0

            print(f"[{cfg.datevshot}] Updating vdrp calibration star maglimit to {magmax}")# and contrast vmin/vmax to {vmin}/{vmax}.")

            files = ['vdrp/vdrp.config','vdrp/vdrp.config.original',
                     'vdrp/vdrp.config.gaia','vdrp/vdrp.config.sdss','vdrp/vdrp.config.panstarrs']

            for file in files:
                system_command(cfg, f"sed -i s/\"mktot_magmax = 22.\"/\"mktot_magmax = {magmax}\"/ {file}")
                #leave the original scaling as is (DD 20260519)
                #system_command(cfg, f"sed -i s/\"cofes_vis_vmin = -20\"/\"cofes_vis_vmin = {vmin}\"/ {file}")
                #system_command(cfg, f"sed -i s/\"cofes_vis_vmax = 30\"/\"cofes_vis_vmax = {vmax}\"/ {file}")

    except:
        print(f"[{cfg.datevshot}] Exception! in update_vdrp_config_limits()", traceback.format_exc())

def node_setup(cfg): #,safelimit=0):
    """

    extra stuff shared for datevshots on same node

    :param cfg:
   # :param safelimit: if positive, do NOT start until the active shot count is BELOW the safelimit
   #                   e.g. THIS shot adds +1 to reach the safelimit
    :return:
    """
    global MaxSafeActiveShots

    try:
        lock = FileLock(Lock_tmp_mutex_fn)
        #if safelimit > 0 and MaxSafeActiveShots > 0:
        if MaxSafeActiveShots > 0:
            redlight = True
            print(f"[{cfg.datevshot}] checking if safe to start ...")

            while redlight:
                with lock:
                    #how many are active?
                    fns = glob.glob(os.path.join(Node_basedir, "*.sync"))
                    active = len(fns)
                    if active < MaxSafeActiveShots: #good to go
                        os.makedirs(Node_basedir, exist_ok=True)
                        with open(os.path.join(Node_basedir, f"{cfg.datevshot}_{os.getpid()}.sync"), "w") as f:
                            #f.write(f"BEGIN {str(datetime.now())}\n")
                            f.write(f"PID: {os.getpid()} BEGAN {str(datetime.now())} from {cfg.cwd_orig}\n")
                        redlight = False
                        print(f"[{cfg.datevshot}] cleared to start.")
                    else:
                        print(f"[{cfg.datevshot}] too many active shots ({active}). Limit {MaxSafeActiveShots}. Must wait ...")

                if redlight:
                    time.sleep(SafeActiveShotsSleep)

        else:
            with lock:
                #create a sync directory and file to hold list of datevshot being used
                #a bit later this will then be the count of simulataneous datevshots and used to tune the num of processes

                os.makedirs(Node_basedir, exist_ok=True)
                with open(os.path.join(Node_basedir,f"{cfg.datevshot}_{os.getpid()}.sync"), "w") as f:
                    f.write(f"PID: {os.getpid()} BEGAN {str(datetime.now())} from {cfg.cwd_orig}\n")
                    #f.write(f"BEGIN {str(datetime.now())}\n")

        # lock auto releases
    except:
        print(f"Exception! in node_setup()",traceback.format_exc())


def node_clean(cfg):
    """

    :param cfg:
    :return:
    """

    if not cfg.node_clean_done:
        try:
            lock = FileLock(Lock_tmp_mutex_fn)
            with lock:
                #clean up my sync
                try:
                    #this might have already been removed (e.g. if post_clean() ran successfully)
                    os.remove(os.path.join(Node_basedir,f"{cfg.datevshot}_{os.getpid()}.sync"))
                    cfg.node_clean_done = True
                except:
                    pass

            # lock auto releases
        except:
            print(f"Exception! in node_clean()",traceback.format_exc())


def node_active_ct(cfg):
    """
    return how many datevshots are active, based on .sync files
    :param cfg:
    :return:
    """

    try:
        ct = -1
        lock = FileLock(Lock_tmp_mutex_fn) #get lock instance, but does not ACQUIRE lock
        with lock: #lock ACQUIRED here
            fns = glob.glob(os.path.join(Node_basedir,"*.sync"))
            ct = len(fns)
            print(f"[{cfg.datevshot}] checking active shot count for shared /tmp: {ct}")
        #lock auto releases
    except:
        print(f"Exception! in node_active_ct()", traceback.format_exc())

    return ct



def num_exposures_in_shot(shotid):
    """

    uses runs<YYYMM> or runt<YYYMM> under Karl's gettar directory

    prior to and including 202407 it is runsYYYYMM
    202408 to 202412 is runs202488
    2025xx is 202500

    :param shotid:
    :return: number of expsoures in the shot and which file was used
    """

    mth = str(shotid)[0:6]
    if int(mth) < 202408:
        fn = os.path.join(karlgettar,f"run?{mth}")
    elif 202408 <= int(mth) <= 202412:
        fn = os.path.join(karlgettar, f"run?202488")
    elif 202500 < int(mth):
        fn = os.path.join(karlgettar, f"run?{mth[0:4]}00")
    else:
        print(f"Invalid shotid? {shotid}")
        return -1, None


    #the ? could be s or t
    try:
        #idx = -1
        date, shot, exp = np.loadtxt(fn.replace('?','s'),usecols=[1,2,3],unpack=True,dtype=str)
        dse = sorted(np.unique([d + s + e for d, s, e in zip(date, shot, exp)]))
        ds = np.array([x[:-5] for x in dse])
        ct = np.count_nonzero(ds == str(shotid))
        if ct > 0:
            fn = fn.replace('?', 's')
    except:
        ct = 0

    #t
    if ct <= 0:
        try:
            ct = 0
            date, shot, exp = np.loadtxt(fn.replace('?', 't'), usecols=[1, 2, 3], unpack=True, dtype=str)
            dse = sorted(np.unique([d + s + e for d, s, e in zip(date, shot, exp)]))
            ds = np.array([x[:-5] for x in dse])
            ct = np.count_nonzero(ds==str(shotid))
            fn = fn.replace('?', 't')
        except:
            pass

    #might not be fatal ... could be a newer observation that is not in the gettars yet, but is still accessible
    # after this call will try to locate the tar file and build a local gettar (make_local_gettar_file)

    #otherwise we did find it ... how may exposures?
    return ct, fn



def get_local_ra_dec(cfg):
    """
    basically, read from the tar file and return the IFU average RA, Dec and track

    :param cfg:
    :return:
    """

    try:
        # since we need other info from the HDU that we do not save, may as well keep
        # it simple and just regenerate the file on a --resume
        # outfile = os.path.join(cfg.cwd, f"{cfg.datevshot}.local_gettar")
        # try:
        #     if os.path.exists(outfile):
        #         print(f"[{cfg.datevshot}] Using existing local gettar file {outfile}")
        #         all_exp = np.loadtxt(outfile,dtype=str,usecols=3,unpack=True)
        #         return len(np.unique(all_exp)) , outfile
        # except:
        #     print(f"[{cfg.datevshot}] Exception using existing local gettar file {outfile}. Will rebuild.",traceback.format_exc())

        #print(f"[{cfg.datevshot}] Building local gettar file for run1s and run2s ... ")

        if cfg.local_het_raw_path is None:
            return None, None, None

        tarfile_path = os.path.join(cfg.local_het_raw_path,
                                    f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")

        tf = tar.open(tarfile_path)
        tarpaths = np.array(tf.getnames())

        sel_fits = np.array([x[-5:] == ".fits" for x in tarpaths])
        tarpaths = tarpaths[sel_fits]

        try:
            u_exp, u_exp_idx = np.unique([x.split("/")[1] for x in tarpaths],return_index=True)
        except: #maybe a bad path? re-run in loop
            l_exp = []
            for x in tarpaths:
                try:
                    l_exp = x.split("/")[1]
                except:
                    print(f"[{cfg.datevshot}] Notice! get_local_ra_dec() unexpected tarpath: {x}")

            u_exp, u_exp_idx = np.unique(l_exp, return_index=True)

        if len(u_exp) == 0: #none found
            print(f"[{cfg.datevshot}] Warning! get_local_ra_dec() could not determine exposures.")
            return None, None, None


        all_ra = []
        all_dec = []
        all_time = []
        track_str = None

        for i in u_exp_idx:
            #just need one for RA, Dec  BUT for time, need one in each exposure
            subtar_fits = tarpaths[i]
            fileh = tf.extractfile(subtar_fits)  # open the fits
            with fits.open(fileh) as hdu:
                # get the cards
                #exp = subtar_fits.split("/")[1]
                #all_exp.append(exp)
                # specid = str(int(hdu[0].header['IFUSLOT'])).zfill(3)
                # ifuslot = str(int(hdu[0].header['SPECID'])).zfill(3)
                # remember this is TRAJRA in decimal hours and we normally want degrees

                all_ra.append(float(hdu[0].header['TRAJRA']) * 15.0) # these should acutally all be the same (also note decimal hours)
                all_dec.append(float(hdu[0].header['TRAJDEC']))

                if float(hdu[0].header['STRUCTAZ']) <= 180.0: #these don't change over exposures
                    cfg.shot_track = 0 #east
                    track_str = 'East'
                else:
                    cfg.shot_track = 1 #west
                    track_str = 'West'

                all_time.append(float(hdu[0].header['PEXPTIME']) - 8.0)

                try: #this does not change per exposures
                    cfg.shot_obj = hdu[0].header['OBJECT'] + " : " + hdu[0].header['QOBJECT'] + \
                               " (" + hdu[0].header['QPROG'] + ")"
                except:
                    pass

            fileh.close()

        if cfg.shot_ra is None or cfg.shot_ra <= -999:
            cfg.shot_ra = np.mean(all_ra)
            cfg.shot_dec = np.mean(all_dec)

        if cfg.total_exp_time is None or cfg.total_exp_time <= 0:
            cfg.total_exp_time = np.sum(all_time)

        print(f"[{cfg.datevshot}] Set shot RA, Dec to ({cfg.shot_ra:0.6f},{cfg.shot_dec:0.6f}) or "
              f"(hours) ({cfg.shot_ra / 15.0:0.6f},{cfg.shot_dec:0.6f}) with {track_str} track")
        print(f"[{cfg.datevshot}] Set shot exp time: {cfg.total_exp_time:0.1f}s")
        print(f"[{cfg.datevshot}] Set shot object: {cfg.shot_obj}")

        tf.close()

    except:
        print(f"[{cfg.datevshot}] Exception! in get_local_ra_dec()", traceback.format_exc())

    return cfg.shot_ra, cfg.shot_dec, cfg.shot_track


def make_local_gettar_file(cfg):
    """
    if we don't have a file from Karl's maveverick/gettar, then just make our own
    using the raw tar file we already copied

     build up a local file that looks like the entries from maverick/gettar/202600 (for example)
     write out that file locally, under this datevshot
     and return the file name
     this will be passed to run1s and run2s

    :param cfg:
    :return: the number of exposures and file name or None
    """

    try:
        #since we need other info from the HDU that we do not save, may as well keep
        # it simple and just regenerate the file on a --resume
        # outfile = os.path.join(cfg.cwd, f"{cfg.datevshot}.local_gettar")
        # try:
        #     if os.path.exists(outfile):
        #         print(f"[{cfg.datevshot}] Using existing local gettar file {outfile}")
        #         all_exp = np.loadtxt(outfile,dtype=str,usecols=3,unpack=True)
        #         return len(np.unique(all_exp)) , outfile
        # except:
        #     print(f"[{cfg.datevshot}] Exception using existing local gettar file {outfile}. Will rebuild.",traceback.format_exc())


        print(f"[{cfg.datevshot}] Building local gettar file for run1s and run2s ... " )

        tarfile_path = os.path.join(cfg.local_het_raw_path,f"{cfg.datevshot[:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar")

        tf = tar.open(tarfile_path)
        tarpaths = np.array(tf.getnames())
        obstype = os.path.basename(tarpaths[0]).split('.')[-2][-3:]
        #only want if is "sci" ... if not, this is an error and the reduction should stop
        if obstype != 'sci':
            print(f"[{cfg.datevshot}] {obstype} != sci")
            tf.close()
            return 0, None

        filestring = ""
        all_exp = []
        all_ra = []
        all_dec = []
        all_track = []
        for subtar_fits in tarpaths:
            fileh = tf.extractfile(subtar_fits)  # open the fits
            with fits.open(fileh) as hdu:
                # get the cards
                exp = subtar_fits.split("/")[1]
                all_exp.append(exp)
                specid = str(int(hdu[0].header['IFUSLOT'])).zfill(3)
                ifuslot = str(int(hdu[0].header['SPECID'])).zfill(3)
                all_ra.append(float(hdu[0].header['TRAJRA'])) #these should acutally all be the same (also note decimal hours)
                all_dec.append(float(hdu[0].header['TRAJDEC']))
                all_track.append(float(hdu[0].header['STRUCTAZ']))

                # example:
                # virus0000009/exp01/virus/20260512T041220.2_013LL_sci.fits',
                # rback 20260131 006 exp01 093 507 202601 0
                # note: the last integer is NOT used in run1s or run2s and is replaced with '2'
                filestring += f"rback {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {exp} {specid} {ifuslot} {str(cfg.datevshot)[0:6]} 2\n"

            fileh.close()

        if len(all_ra) > 0:
            cfg.shot_ra = np.nanmean(all_ra) * 15.0 #remember this is TRAJRA in decimal hours and we normally want degrees
            cfg.shot_dec = np.nanmean(all_dec)
            az = np.nanmean(all_track)
            if az <= 180.0:
                cfg.shot_track = 0 #east
            else:
                cfg.shot_track = 1  #west
            print(f"[{cfg.datevshot}] Set shot RA, Dec to ({cfg.shot_ra:0.6f},{cfg.shot_dec:0.6f}) or (hours) ({cfg.shot_ra/15.0:0.6f},{cfg.shot_dec:0.6f})")

        tf.close()

        if len(filestring) > 1:
            outfile = os.path.join(cfg.cwd,f"{cfg.datevshot}.local_gettar")
            with open(outfile, "w") as fgt:
                for fs in filestring:
                    fgt.write(fs)

            print(f"[{cfg.datevshot}] Wrote {outfile} for run1s and run2s ... ")
            cfg.gettar_fn = outfile
            return len(np.unique(all_exp)), outfile

    except:
        print(f"[{cfg.datevshot}] Exception! in make_local_gettar_file()", traceback.format_exc())
        return 0, None

def run_run1s(cfg):
    """

    science_reduction/sciscripts/run1  batch file


    :return:
    """

    print(f"[{cfg.datevshot}] run1s ...")
    #probably should change to use subprocess
    if cfg.exp > 0:
        exps = [cfg.exp]
    else:
        exps = np.arange(1,cfg.numexp+1)

    for exp in exps:
        #example run1s 20240730 009 exp01 202407

        system_command(cfg,f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# run1s") #use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# run2s")  # use '#' as sed separator rather than "/"
        #system_command(cfg,f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# rtaremc")  # use '#' as sed separator rather than "/"

        if cfg.multifits_only:
            if cfg.sub_shot is not None:
                if cfg.ifuslot is not None:
                    print(f"[{cfg.datevshot}] substituting shot {cfg.sub_shot} for gettar lookup and using only IFUSlot {cfg.ifuslot}...")
                    system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 {cfg.sub_shot} $3 {cfg.ifuslot}'# run1s")
                    system_command(cfg, f"sed -i s#'$1,$2,$3,$4'#'$1,$2,\"{cfg.datevshot[-3:]}\",$4'# run1s")
                    system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 {cfg.sub_shot} $3 {cfg.ifuslot}'# run2s")
                    system_command(cfg, f"sed -i s#'$1,$2,$3,$4'#'$1,$2,\"{cfg.datevshot[-3:]}\",$4'# run2s")
                else:
                    print(f"[{cfg.datevshot}] substituting shot {cfg.sub_shot} for gettar lookup ...")
                    system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 {cfg.sub_shot} $3'# run1s")
                    system_command(cfg, f"sed -i s#'$1,$2,$3,$4'#'$1,$2,\"{cfg.datevshot[-3:]}\",$4'# run1s")
                    system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 {cfg.sub_shot} $3'# run2s")
                    system_command(cfg, f"sed -i s#'$1,$2,$3,$4'#'$1,$2,\"{cfg.datevshot[-3:]}\",$4'# run2s")
            elif cfg.ifuslot is not None:
                print(f"[{cfg.datevshot}] using only IFUSlot {cfg.ifuslot}...")
                system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 $2 $3 {cfg.ifuslot}'# run1s")
                system_command(cfg, f"sed -i s#'$1 $2 $3'#'$1 $2 $3 {cfg.ifuslot}'# run2s")

        #cmd = "sed -i s#\${scriptdir}"+f"#{cfg.scriptdir}/sciscripts/# run1s"
        #scripts have already been copied to shot workding dir
        cmd = "sed -i s#\${scriptdir}" + f"#{cfg.cwd}/# run1s"
        system_command(cfg,cmd)  # use '#' as sed separator rather than "/"

        # awk '{print $1,$2,$3,$4,$5,$6,$7,2}' > rj
        #make all the calls background so runs the whole set if IFUs at once
        # match_str = "2}' > rj"
        # sub_str ="2,\"\&\"}' > rj"
        # system_command(cfg, f"sed -i s/\"{match_str}\"/\"{sub_str}\"/ run1s")

        #actually run it here
        system_command(cfg,f"run1s {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} exp{str(exp).zfill(2)} {cfg.datevshot[0:6]}")


    #need this later for vdrp pathing
    system_command(cfg,"ln -s ../reductions reductions")


def check_run1s(cfg):
    """
    todo: any automated checks that can be performed
          correct any possible
          log non-fatal issues
          terminate if fatal

          Previously this is handled in a notebook on my machine: calibrations/check_red1.ipynb


    :param cfg:
    :return:
    """

    print(f"[{cfg.datevshot}] todo:  check run1s ... see check_red1.ipynb")

    rc = 0

    return rc


def vdrp_check_norms(cfg):
    """
    basic check on dither norms (if there are any
    :param cfg:
    :return:
    """

    print("check norms ... ")
    #only run IF there are dithers (multiple exposures ... assume dither)
    if cfg.exp == 0 and cfg.numexp > 1:
        rc = 0
        vdrp_path = "./"

        try:
            fns = glob.glob(os.path.join(vdrp_path, f"{cfg.datevshot[0:6]}*v???"))
        except:
            fns = glob.glob(os.path.join(vdrp_path, "20*v???"))

        if len(fns) == 0:
            rc = -1
            print(f"[{cfg.datevshot}] vdrp : check norms ... no matches found")
        else:
            fns = sorted(fns)
            for fn in fns:
                norms = np.loadtxt(os.path.join(fn, "norm.dat"))  # one line, 3 values
                cfg.dither_norms = norms
                if np.count_nonzero(abs(1 - norms) > 0.5) > 0 or np.any(np.isnan(norms)):
                    print("Possible bad dither norm:", os.path.basename(fn), norms)
                    rc = -1

        if rc < 0:
            print("check norms ... fail")
        else:
            print("check norms ... pass")
        return rc
    else:
        print("check norms ... OK")
        return 0

def vdrp_check_shout_ifu(cfg):
    """

    :param cfg:
    :return:
    """


    # usually a YYYYMM
   # cl_args = list(map(str.lower, sys.argv))

    path = "./"  # /scratch/03261/polonius/science_reductions/vdrp/shifts/"
    rc = 0

    print(f"[{cfg.datevshot}] vdrp : check shout ifu ... ")
    try:
        wildcard = cfg.datevshot[0:6]
    except:
        wildcard = ""

    paths = glob.glob(f"{path}{wildcard}*v???")

    print(f"[{cfg.datevshot}] vdrp : checking {len(paths)} directories ... ")
    if len(paths) == 0:
        rc = -1
        print(f"[{cfg.datevshot}] vdrp : check shout ifu ... no matches found for: {path}{wildcard}*v??? ")

    # for path in tqdm(paths):
    for path in paths:
        basedir = None
        try:
            fn = os.path.join(path, "shout.ifu")

            basedir = os.path.abspath(fn).split("/")[-2]

            stats = os.stat(fn)

            if stats.st_size == 0:
                print(f"[{cfg.datevshot}] vdrp : {basedir} : empty shout.ifu")
                rc = -1
            elif stats.st_size < 1000:
                print(f"[{cfg.datevshot}] vdrp : {basedir} : small shout.ifu ({stats.st_size})")
                rc = -1

        except:
            print(traceback.format_exc())
            print(f"[{cfg.datevshot}] cvdrp : {basedir} : unknown or missing shout.ifu")

    if rc != 0:
        print(f"[{cfg.datevshot}] vdrp : check shout ifu ... fail ")
    else:
        print(f"[{cfg.datevshot}] vdrp : check shout ifu ... OK")
    return rc

def vdrp_cp2dithall(cfg,catalog=None):
    """

    :param cfg:
    :return:
    """

    print(f"[{cfg.datevshot}] vdrp : cp2dithall ... ")
    rc = 0
    dirs = glob.glob(f"./{cfg.datevshot[0:6]}??v???")
    #note: do not want the log file YYYYMMDDvSSS.log
    if catalog is None:
        dithall_dir = "dithall"
    else:
        dithall_dir = "dithall."+catalog

    for d in dirs:
        try:
            dithall_use = os.path.join(d, "dithall.use")
            if not os.path.exists(dithall_use):
                rc = -1
                print(f"[{cfg.datevshot}] vdrp : cp2dithall {catalog} ... fail. File does not exist {os.getcwd()} {dithall_use}")
                continue

            with open(dithall_use, "r") as f1:
                # skip 1st line
                _ = f1.readline()
                os.makedirs(dithall_dir,exist_ok=True)
                outfn = os.path.basename(d)
                outfn = os.path.join("./", f"{dithall_dir}/{outfn}.dithall")
                with open(outfn, "w+") as f2:
                    for line in f1:
                        f2.write(line)
        except:
            print(traceback.format_exc())
            rc = -1

    if rc != 0:
        print(f"[{cfg.datevshot}] vdrp : cp2dithall {catalog}... fail")
    else:
        print(f"[{cfg.datevshot}] vdrp : cp2dithall {catalog}... OK")
    return rc

def run_vdrp(cfg):
    """

    :param cfg:
    :return:
    """

    def do_gaia(cfg):
        #command is based on rta.YYYYMM
        #run_shifts 20240730 009 16.317927 33.689304 1  20240730 GAIA
        fail_gaia = False
        print(f"[{cfg.datevshot}] VDRP: GAIA")
        try:
            with open(f"rta.{cfg.datevshot[0:6]}", "r") as rta:
                for line in rta:  # really should only be one line
                    if len(line) > 10:
                        # we DON'T want the carriage return
                        line = line.rstrip('\n')
                        system_command(cfg, f"{line} {cfg.datevshot[0:8]} GAIA")

            vdrp_check_norms(cfg)
            vdrp_check_shout_ifu(cfg)
            vdrp_cp2dithall(cfg, "gaia")
            # move the output under gaia
            os.makedirs("gaia", exist_ok=True)
            system_command(cfg, f"mv {cfg.datevshot[0:6]}??v??? gaia")

            # check the dithall.gaia exists
            if not os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.gaia")):
                fail_gaia = True

        except Exception as e:
            print(f"[{cfg.datevshot}] VDRP: GAIA fail.", e, "\n", traceback.format_exc())
            fail_gaia = True

        return fail_gaia

    def do_sdss(cfg):
        #command is based on rta.YYYYMM
        #run_shifts 20240730 009 16.317927 33.689304 1  20240730 SDSS
        fail_sdss = False
        print(f"[{cfg.datevshot}] VDRP: SDSS")
        try:
            with open(f"rta.{cfg.datevshot[0:6]}", "r") as rta:
                for line in rta:  # really should only be one line
                    if len(line) > 10:
                        # we DON'T want the carriage return
                        line = line.rstrip('\n')
                        system_command(cfg, f"{line} {cfg.datevshot[0:8]} SDSS")

            vdrp_check_norms(cfg)
            vdrp_check_shout_ifu(cfg)
            vdrp_cp2dithall(cfg, "sdss")
            # move the output under sdss
            os.makedirs("sdss", exist_ok=True)
            system_command(cfg, f"mv {cfg.datevshot[0:6]}??v??? sdss")

            # check the dithall.sdss exists
            if not os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.sdss")):
                fail_sdss = True

        except Exception as e:
            print(f"[{cfg.datevshot}] VDRP: SDSS fail.", e, "\n", traceback.format_exc())
            fail_sdss = True

        return fail_sdss

    #change to the vdrp/shifts directory
    # since we do not know if we will need all of GAIA, SDSS and PanSTARRS at this point,
    #    run all three of them


    os.chdir(os.path.join(cfg.cwd,"vdrp/shifts"))

    rta_fn = f"rta.{cfg.datevshot[0:6]}"

    if cfg.build_rta:
        #track is east (0) or west (1)
        with open(rta_fn, "w") as f:
            f.write(f"run_shifts.sh {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {cfg.shot_ra / 15.0} {cfg.shot_dec} {cfg.shot_track} \n")


    if not os.path.exists(rta_fn):
        if os.path.exists(os.path.join(karlgettar,rta_fn)):
            shutil.copy(os.path.join(karlgettar,rta_fn), ".")

    #one last try: either the copy failed or the auto-creation failed
    if not os.path.exists(rta_fn):
        with open(rta_fn, "w") as f:
            f.write(f"run_shifts.sh {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {cfg.shot_ra / 15.0} {cfg.shot_dec} {cfg.shot_track} \n")

    #if not os.path.exists(rta_fn): #this is fatal now, but it will fail appropriately

    fail_gaia = False
    fail_sdss = False
    fail_panstarrs = False
    #not quite the opposite of the fail_xxx
    used_gaia = False
    used_sdss = False
    used_panstarrs = False


    if cfg.starcat_ast == None or cfg.starcat_ast == 'gaia':
        #GAIA first
        fail_gaia = do_gaia(cfg)
        if fail_gaia:
            print(f"[{cfg.datevshot}] VDRP: Astrometry FAIL with GAIA. Attempting SDSS as fallback.")
            fail_sdss = do_sdss(cfg)
            if fail_sdss:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback FAIL with SDSS.")
            else:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback PASS with SDSS.")
                used_sdss = True
        else:
            used_gaia = True
    elif cfg.starcat_ast == 'sdss':
        fail_sdss = do_sdss(cfg)
        if fail_sdss:
            print(f"[{cfg.datevshot}] VDRP: Astrometry FAIL with SDSS. Attempting GAIA as fallback.")
            fail_gaia = do_sdss(cfg)
            if fail_gaia:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback FAIL with GAIA.")
            else:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback PASS with GAIA.")
                used_gaia = True
        else:
            used_sdss = True
    else: #should not happend
        print(f"[{cfg.datevshot}] VDRP: unexpected starcat_ast value [{cfg.starcat_ast}] ; using GAIA as default ")
        fail_gaia = do_gaia(cfg)
        if fail_gaia:
            print(f"[{cfg.datevshot}] VDRP: Astrometry FAIL with GAIA (2). Attempting SDSS as fallback.")
            fail_sdss = do_sdss(cfg)
            if fail_sdss:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback FAIL with SDSS (2).")
            else:
                print(f"[{cfg.datevshot}] VDRP: Astrometry fallback PASS with SDSS (2).")
                used_sdss = True
        else:
            used_gaia = True


    #!!! notice: we are doing limited checking here ... that will be done later
    # previously would run: check_norms YYYYMM   (included)
    #                       check_shot.ifu YYYYMM (included)
    #  examine the .pngs manually  NOT DONE
    #                  run: make_good_shots YYYYMM  NOT DONE
    #                       cp2dithall YYYYMM  (included)


    #always third
    if do_panstarrs or (fail_gaia and fail_sdss):
        #PanSTARRS
        print(f"[{cfg.datevshot}] VDRP: PANSTARRS")
        try:
            with open(f"rta.{cfg.datevshot[0:6]}","r") as rta:
                for line in rta: #really should only be one line
                    if len(line) > 10:
                        #we DON'T want the carriage return
                        line = line.rstrip('\n')
                        system_command(cfg,f"{line} {cfg.datevshot[0:8]} PANSTARRS")

            vdrp_check_norms(cfg)
            vdrp_check_shout_ifu(cfg)
            vdrp_cp2dithall(cfg,"panstarrs")
            #move the output under panstarrs
            os.makedirs("panstarrs",exist_ok=True)
            system_command(cfg,f"mv {cfg.datevshot[0:6]}??v??? panstarrs")

            # check the dithall.panstarrs exists
            if not os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.panstarrs")):
                fail_panstarrs = True
            else:
                used_panstarrs = True

        except Exception as e:
            print(f"[{cfg.datevshot}] VDRP: PANSTARRS fail.", e, "\n", traceback.format_exc())
            fail_panstarrs = True

    #update to which one was actually used (should be 0 or 1 at most)
    if used_gaia:
        if cfg.starcat_ast != 'gaia':
            print(f"[{cfg.datevshot}] VDRP: updating astrometry catalog used to GAIA")
            cfg.starcat_ast = "gaia"
    elif used_sdss:
        if cfg.starcat_ast != 'sdss':
            print(f"[{cfg.datevshot}] VDRP: updating astrometry catalog used to SDSS")
            cfg.starcat_ast = "sdss"
    elif used_panstarrs:
        if cfg.starcat_ast != 'panstarrs':
            print(f"[{cfg.datevshot}] VDRP: updating astrometry catalog used to PANSTARRS")
            cfg.starcat_ast = "panstarrs"
    else: #going to be fatal
        print(f"[{cfg.datevshot}] VDRP: unable to fix astrometry.")

    #todo: which is the main dithall, etc???? (normally it is SDSS for calibration and gaia for astrometry)
    #print("!!! todo: copy GAIA dithaall to /scatch/projects and /corral-repl ???")

    # set the star catalog to use for calibration (make a softlink to the starcat specific output for this shot)
    # note: this happens regardless of outcome
    # just in case something went badly wrong, make sure we are in the right directory
    os.chdir(os.path.join(cfg.cwd, "vdrp/shifts"))
    if os.path.exists(f"{cfg.starcat_cal}/{cfg.datevshot}"):
        system_command(cfg, f"ln -s {cfg.starcat_cal}/{cfg.datevshot} {cfg.datevshot}")

    os.chdir(cfg.cwd)


def check_vdrp(cfg):
    """

    :param cfg:
    :return:
    """

    fail_gaia = False
    fail_sdss = False
    fail_panstarrs = False

    # check the dithall.gaia exists
    if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.gaia")):
        print(f"[{cfg.datevshot}] VDRP: GAIA [Pass]")
    else:
        print(f"[{cfg.datevshot}] VDRP: GAIA [FAIL]")
        fail_gaia = True

    # check the dithall.sdss exists
    if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.sdss")):
        print(f"[{cfg.datevshot}] VDRP: SDSS [Pass]")
    else:
        print(f"[{cfg.datevshot}] VDRP: SDSS [FAIL]")
        fail_sdss = True

    # check the dithall.panstarrs exists (this might not have been run, so it may be okay if it did not)
    if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.panstarrs")):
        print(f"[{cfg.datevshot}] VDRP: PanSTARRs [Pass]")
    else:
        fail_panstarrs = True
        if (fail_gaia and fail_sdss):
            print(f"[{cfg.datevshot}] VDRP: PanSTARRs [FAIL]")
        else:
            #maybe it just was not set to run
            print(f"[{cfg.datevshot}] VDRP: PanSTARRs not found.")

    if fail_gaia and fail_sdss and fail_panstarrs:
        return -1
    else:
        return 0



def prepare_reduction_dir(cfg):
    """
    basically create the expected paths and copy the reduction (multi*fits) files
    where they are expected

    :param cfg:
    :return:
    """

    try:
        rc = 0

        os.chdir(cfg.cwd)
        expdirs = glob.glob("d20*exp??")

        #/scratch/03261/polonius/red1/reductions/20241021/virus/virus0000007/exp01/virus/multi_501_080_012_RL.fits

        for expdir in expdirs:
            exp = os.path.basename(expdir)[-5:]
            virus_shot = "virus0000" + os.path.basename(expdir)[-8:-5]
            date = cfg.datevshot[0:8]
            if red1path is not None:
                datadir = os.path.join(red1path,f"{date}/virus/{virus_shot}/{exp}/virus")
            else:
                datadir = os.path.join(cfg.cwd_orig, f"reductions/{date}/virus/{virus_shot}/{exp}/virus")

            tarfile = f"d{cfg.datevshot[0:8]}s{cfg.datevshot[-3:]}{exp}_mu.tar"

            print(f"[{cfg.datevshot}] Creating directory and untarring ({expdir}/{tarfile}) multi*fits to: {datadir}")

            Path(datadir).mkdir(parents=True, exist_ok=True)

            #mutlti*.fits
            cmd = f"tar -xvf {expdir}/{tarfile} -C {datadir}"
            system_command(cfg,cmd)

            #CoFe*.fits
            tarfile = f"d{cfg.datevshot[0:8]}s{cfg.datevshot[-3:]}{exp}_co.tar"
            print(f"[{cfg.datevshot}] Untarring ({expdir}/{tarfile}) CoFe*fits to: {datadir}")
            cmd = f"tar -xvf {expdir}/{tarfile} -C {datadir}"
            system_command(cfg,cmd)

    except:
        print(traceback.format_exc())
        rc = -1
        print(f"Fatal. Could not prepare reductions directory.")

    return rc



def run_fluxcalibration(cfg,star_catalog='sdss'):
    """

    :param cfg:
    :return:
    """

    # NOTICE: this operates under two directories ... there is a second directory change partway down

    print(f"[{cfg.datevshot}] flux calibration using {star_catalog} ... ")

    #setup the softlink for the star_catalog
    os.chdir(os.path.join(cfg.cwd, "vdrp/shifts"))
    if os.path.exists(cfg.datevshot):
        system_command(cfg, f"unlink {cfg.datevshot}")
    if os.path.exists("dithall"):
        system_command(cfg, "unlink dithall")

    system_command(cfg, f"ln -s {star_catalog}/{cfg.datevshot} {cfg.datevshot}")
    system_command(cfg, f"ln -s dithall.{star_catalog} dithall")

    os.chdir(os.path.join(cfg.cwd,"detect"))

    #copy dithal to detect dir so can be softlinked in rsp3fc
    #this is for use under Karl's rsetstar ... in rsp3fc a softlink is created for this .dithall file
    #under the /tmp/rsXXXX working dir for rsetstar (really ~gebhardt/bin/fitradecsp)
    #s|t if the usual .dithall for the datevshot under /scratch/projects/hetdex/detect is NOT found (e.g. we are post-HETDEX)
    #then it will default to the local.dithall
    if len(glob.glob("*.dithall")) > 0:
        system_command(cfg, f"unlink *.dithall")
    system_command(cfg, f"ln -s {os.path.join(cfg.cwd, 'vdrp/shifts/dithall/*.dithall')} .")

    if os.path.exists("local.tp"):
        system_command(cfg, f"unlink local.tp")
    system_command(cfg,f"ln -s tp/{cfg.datevshot}sedtp_f.dat local.tp")

    print(f"[{cfg.datevshot}] flux calibration using {star_catalog} (rsetstar) ... ")
    #call rsetstar independently
    system_command(cfg, f"rsetstar {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {star_catalog}")

    print(f"[{cfg.datevshot}] flux calibration using {star_catalog} (rallcal) ... ")


    #may need to update the tp upper limit
    #for HETDEX shots that are long, I think this is still okay, since it just alters the upper limit
    #and those shots that go longer are usually due to poor seeing, which pushes the other way?


    if cfg.total_exp_time is not None:
        # numexp = max(1,cfg.numexp)
        # mux = max(1.0, (cfg.total_exp_time / numexp / 360.)) #360secs is the nominal default for a HETDEX exposure
        # tp_upper = max(0.25, 0.25 * mux)

        tp_upper = max(0.25, 0.25 * linear_exptime_scale(cfg))

        if tp_upper > 0.25:
            print(f"[{cfg.datevshot}] *** Altering max tp from 0.25 to {round(tp_upper,2)} due to long exptime ({cfg.total_exp_time}s)")
            base_tp_str = "$2<0.25"
            tp_str = f"$2<{round(tp_upper,2)}"

            #we are under detect currently
            system_command(cfg, f"sed -i s/\"{base_tp_str}\"/\"{tp_str}\"/ rgettp")
            system_command(cfg, f"sed -i s/\"{base_tp_str}\"/\"{tp_str}\"/ cal_script/rgettp")
        else:
            print(f"[{cfg.datevshot}] Exposure is typical. Will not adjust tp calculation.")
    else:
        print(f"[{cfg.datevshot}] *** WARNING!!! exposure time is not known. Cannot adjust tp calculation.")

    #no longer includes rsetstar
    system_command(cfg,f"rallcal {cfg.datevshot[0:8]} {cfg.datevshot[-3:]}")

    return 0


def check_fluxcalibration(cfg):
    """

    :param cfg:
    :return:
    """
    #print("todo: check flux calibration ... ")
    rc = 0
    os.chdir(os.path.join(cfg.cwd, "detect"))

    #todo: if tp/*setp_.dat is "bad" (last columns all zero)
    #      the try with a different star catalog ... e.g.
    #      by default, SDSS is used for calibration, but if that fails
    #      try GAIA and if that fails
    #      try PANSTARRS
    #          >> to do that, need to change the softlink under vdrp/shifts/YYYYMMDDsSSS to gaia/YYYYMMDDvSSS or panstarrs/YYYYMMDDvSSS
    #          >> if they exist (e.g. that might have failed during vdrp)
    #          >> then, rerun   run_fluxcalibration() and check again

    out = None
    try:
        #20240730v006sedtp_f.dat
        if os.path.exists(f"tp/{cfg.datevshot}sedtp_f.dat"):
            out = np.loadtxt(f"tp/{cfg.datevshot}sedtp_f.dat")
            if np.shape(out)[0] > 0:
                if not np.any(out[:, 5]):  # 5 is the actual throughput, I think and 4 is tied to it? cols 2, 3 don't seem to matter
                    rc = -1
                    print(f"[{cfg.datevshot}] bad throughput {cfg.datevshot}. All zero values.")
            else:
                rc = -1
                print(f"[{cfg.datevshot}] bad throughput {cfg.datevshot}. *sedtp_f.dat is empty.")
        else:
            rc = -1
            print(f"[{cfg.datevshot}] bad throughput {cfg.datevshot}. *sedtp_f.dat file does not exist.")
    except:
        print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] bad throughput {cfg.datevshot}")
        print(f"[{cfg.datevshot}] bad throughput. output shape = {np.shape(out)}")

    return rc

#############################
#step 04
#############################




def run_make_ifucen(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        os.chdir(os.path.join(cfg.cwd, "getcen"))
        system_command(cfg, f"rgetifucen {cfg.datevshot[0:8]} {cfg.datevshot[-3:]}")

        #not much to check here
        if os.path.exists(f"ifucen_{cfg.datevshot}.dat"):

            #check the RA, decs
            ras, decs = np.loadtxt(f"ifucen_{cfg.datevshot}.dat", dtype=float,
                                   usecols=(1, 2), unpack=True)
            badct = np.count_nonzero(ras == -666)
            if badct > 0:
                if badct < 3:  # allow 1 or 2 bad ifus
                    badsel = ras == -666
                    mf = np.loadtxt(f"ifucen_{cfg.datevshot}.dat", dtype=str, usecols=(0), unpack=True)
                    print(f"[{cfg.datevshot}] failed to get {badct} IFU centers ({mf[badsel]}). Treat as non-fatal.")
                    rc = 0 #call it non-fatal
                else:
                    print(f"[{cfg.datevshot}] failed to get {badct} IFU centers. Likely fatal.")
                    rc = -1
            else:
                #all is good, assume anyway
                rc = 0
        else:
            rc = -1

    except:
        print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] failed to get IFU centers.")


    return rc




def run_rfft(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        rc = 0
        os.chdir(os.path.join(cfg.cwd, "alldet"))

        output_suffixes = ["amp.dat","ds9.reg","sky.dat","sub.fits"]

        #iterate over exp
        if cfg.exp > 0:
            exps = [cfg.exp]
        else:
            exps = np.arange(1, cfg.numexp + 1)

        for exp in exps:
            system_command(cfg, f"rfft {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} exp{str(exp).zfill(2)}")

            # check the output direcotry for 4 files per exp
            # *amp.dat, *ds9.reg, *sky.dat, *sub.fits
            # like: d20230904s011exp01amp.dat
            for suffix in output_suffixes:
                fn = f"output/d{cfg.datevshot.replace('v','s')}exp{str(exp).zfill(2)}{suffix}"
                if os.path.exists(fn):
                    if suffix == "amp.dat":
                        #todo: check certain columns
                        col1 = np.loadtxt(fn, dtype=str, usecols=[0],skiprows=1)
                        #known bad pattern
                        #IF you see this, it may be a problem with the "list" of multi*fits for rfft
                        #  (see alldet/rfft and the creation of "list") ... there is a 120 character limit
                        #  in the Fortran code for each list entry and previously this issue occured when
                        #  the entries were 166 characters each
                        bad_pattern = '\x00\x00\x00_\x00\x00\x00_\x00\x00\x00_'
                        if col1[0] == bad_pattern:
                            print(f"[{cfg.datevshot}] rfft [FAIL] {fn}. Bad data/format.")
                        else:
                            #assume good
                            print(f"[{cfg.datevshot}] rfft [Pass] {fn}")
                    elif suffix == "ds9.reg": #not important, but just assume good
                        print(f"[{cfg.datevshot}] rfft [Pass] {fn}")
                    elif suffix == "sky.dat": #this is the easy one to check
                        col2 = np.loadtxt(fn,dtype=float,usecols=[1])
                        if np.any(col2):
                            print(f"[{cfg.datevshot}] rfft [Pass] {fn}")
                        else:
                            print(f"[{cfg.datevshot}] rfft [FAIL] {fn}. All zero values.")
                    else: #this is the fits ... do we want to check it?
                        print(f"[{cfg.datevshot}] rfft [Pass] {fn}")

                else:
                    print(f"[{cfg.datevshot}] rfft [FAIL] {fn} file not found")
                    rc = -1 #even though this is fatal, go ahead and loop over all files and expXX so can get into the log
                            #could be a useful diagnoistic

    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Exception! run_rfft", traceback.format_exc())

    avg_sky = get_avg_sky(cfg)
    if avg_sky is not None:
        cfg.avg_sky = avg_sky
        if avg_sky > MAX_SAFE_AVG_SKY:
            if avg_sky >= FAIL_AVG_SKY:
                print(f"[{cfg.datevshot}] Average Sky is catastrophically large and/or unable to be fit: {avg_sky:0.1f}.")
                rc = -1  # fatal
            else:
                #not fatal, but still a warning
                print(f"[{cfg.datevshot}] Average Sky is problematically large: {avg_sky:0.1f}. Non-fatal.")
    return rc



def mp_rcal_worker(out_list,cfg,set_idx,indicies,multis, ras, decs):
    """
    worker for multithreaded calls to run_rcal

    :param cfg:
    :param out_list: if tracking results, they would go here
    :return:
    """

    print(f"[{cfg.datevshot}] idx[{set_idx}] mp_rcal_worker serial for: {multis}",flush=True)
    for multi, ra, dec, ix in zip(multis, ras, decs,indicies):
        print(f"[{cfg.datevshot}] (run_rcal) {set_idx}-{ix}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ",flush=True)  # ,end="")

        # grid_n = 3  # n x n so 3 = a 3x3 grid
        # grid_step = 0.5  # grid step size
        #
        # if cfg.numexp < 3 or cfg.dither_configuration == 0:
        #     grid_n = 7
        #NOTE: the 0.5 and 3 here are not actually used, so do not update
        cmd = f"rcal_all {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 106"
        system_command(cfg, cmd)


def mp_rcal(cfg,multis, ras, decs,num_procs=NumProcs_mp_rcal):
    """

    these are running basically blind ... the output is to files and I only need to know that it is done,
    regardless of success or fail

    so, break up in to XXX tasks here and spin up, then return when all done

    :param cfg:
    :return:
    """

    #print(f"Top of mp_rcal: len(multis) = {len(multis)}")

    #tune num procs
    active_shots = node_active_ct(cfg)
    if active_shots > 0:
        num_procs = int(np.floor(MaxTotalProcs_mp_rcal / active_shots))
        num_procs = max(1,min(num_procs,MaxPerShotProcs_mp_rcal)) #always at least one

    print(f"[{cfg.datevshot}] flux calibration  ***MULTITHREADED*** (run_rcal) ({len(ras)}), num_procs = {num_procs}...")

    with multiprocessing.Manager() as manager:
        mgr_list = manager.list()
        idx = np.array_split(np.arange(len(multis)),num_procs)
        processes = []

        for i in range(num_procs):
            print(f"[{cfg.datevshot}] (run_rcal): Spinning up mp_rcal i={i},  len(multis) = {len(multis[idx[i]])}, idx[i] = {idx[i]}")
            process = multiprocessing.Process(target=mp_rcal_worker,
                                               args=(mgr_list,cfg,i,idx[i],multis[idx[i]],ras[idx[i]],decs[idx[i]]))

            processes.append(process)
            process.start()

        # Wait for processes to complete
        for process in processes:
            print(f"[{cfg.datevshot}] (run_rcal): Joining {process}")
            process.join()

        print(f"[{cfg.datevshot}] (run_rcal): mp_rcal all done")
        #again, don't care about the results here, just need them done

def run_rcal(cfg):
    """

    :param cfg:
    :return:
    """

    #currently bullet 3,4 under document

    #rcal_all 244.7383880 33.8140602 35 4505 50 512_104_026 20240730v009 1.70 3.0 3.5 0.5 3 106
    #rcal_all 244.7716830 33.8165741 35 4505 50 513_105_051 20240730v009 1.70 3.0 3.5 0.5 3 106
    #rcal_all 244.7050930 33.8114395 35 4505 50 514_103_019 20240730v009 1.70 3.0 3.5 0.5 3 106

    try:
        failed_rcal_list = []
        passed_rcal_list = []
        rc = 0
        os.chdir(os.path.join(cfg.cwd, "alldet")) #make sure we are in the right directory
        os.makedirs("cal_out",exist_ok=True)


        #prep, needed for rcal_all copy to /tmp
        system_command(cfg, f"cp ../detect/{cfg.datevshot}/norm.dat .")
        system_command(cfg, f"cp ../detect/{cfg.datevshot}/fwhm.out fwhm.use")

        multis = np.loadtxt(os.path.join("../getcen/",f"ifucen_{cfg.datevshot}.dat"),dtype=str,usecols=(0),unpack=True)
        ras,decs = np.loadtxt(os.path.join("../getcen/", f"ifucen_{cfg.datevshot}.dat"), dtype=float,usecols=(1,2),unpack=True)
        # print(f"run_rcal ({len(ras)})...")
        # ct = 0
        # for multi,ra,dec in zip(multis,ras,decs):
        #     ct +=1
        #     print(f"{ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ",end="")
        #     system_command(cfg,f"rcal_all {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 106")

        #print(f"run_rcal ***MULTITHREADED*** ({len(ras)})...")
        mp_rcal(cfg, multis, ras, decs)#, num_procs=6)

        #now check them all
        ct = 0
        for multi, ra, dec in zip(multis, ras, decs):
            ct += 1
            #check the output exists cal_out/20240730v009_514_103_019_cal.fits
            base_str = f"[{cfg.datevshot}] checking run_cal ({ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... "
            #print(f"{ct}) checking {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ", end="")
            outfn = f"{cfg.datevshot}_{multi[6:]}_cal.fits"
            if os.path.exists(os.path.join("cal_out/", outfn)):
                #check the contents
                try:
                    hdu = fits.open(os.path.join("cal_out/", outfn))
                    #the 5 HUDs:  [0] calib, [1] calibe, [2] calib_c, [3] calibe_c, [4] fullsky
                    #data is (usuall) 1344x1036 ... 112 fibers x 4 amps x 3 dithers by 1036 rectified wavelength bins
                    if len(hdu) != 5:
                        print(f"{base_str} [FAIL]: bad hdu length. got {len(hdu)}, expected 5")
                        failed_rcal_list.append(outfn)
                    else:
                        all_zero_list = []
                        for i in range(5):
                            if np.any(hdu[i].data):
                                pass #assume otherwise okay
                            else:
                                all_zero_list.append(i)

                        if len(all_zero_list) > 0:
                            print(f"{base_str} [FAIL]: all zeroes for HDU{all_zero_list}")
                            failed_rcal_list.append(outfn)
                        else:
                            print(f"{base_str}  [Pass]")
                            passed_rcal_list.append(outfn)
                except:
                    print(traceback.format_exc())
                    print(f"{base_str}  [FAIL]. Exception!")
                    failed_rcal_list.append(outfn)

            else:
                print(f"{base_str}  [FAIL]: no cal_out")
                failed_rcal_list.append(outfn)

        if len(passed_rcal_list) == len(ras):
            rc = 0
            print(f"[{cfg.datevshot}] (run_rcal) All Pass")
        elif len(failed_rcal_list) == len(ras):  # all failed
            rc = -1
            print(f"[{cfg.datevshot}] (run_rcal) ALL FAIL")
        else:
            rc = 1
            print(f"[{cfg.datevshot}] (run_rcal) Mixed results of {len(ras)}: {len(passed_rcal_list)} Pass, {len(failed_rcal_list)} FAIL")


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"{cfg.datevshot} (run_rcal): Fatal exception!",traceback.format_exc())

    return rc


def mp_rf1_worker(out_list,cfg,set_idx,indicies,multis, ras, decs):
    """
    worker for multithreaded calls to rf1

    :param cfg:
    :param out_list: if tracking results, they would go here
    :return:
    """

    print(f"[{cfg.datevshot}] idx[{set_idx}] mp_rf1_worker serial for: {multis}",flush=True)
    for multi, ra, dec, ix in zip(multis, ras, decs,indicies):
        print(f"[{cfg.datevshot}] (mp_rf1) : {set_idx}-{ix}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ",flush=True)  # ,end="")

        #grid_n = 3 #n x n so 3 = a 3x3 grid
        #grid_step = 0.5 #grid step size
        #print(f"[{cfg.datevshot}] *** test *** set grid_step and grid_n to 0.0 1")
        # grid_n = 1
        # grid_step = 0.0
        grid_n = f"{int(cfg.linedet_parms[0])}"
        grid_step = f"{float(cfg.linedet_parms[1]):0.1f}"
        #if cfg.numexp < 3 or cfg.dither_configuration == 0:
        #    grid_n = 13
        #    grid_step = 0.25

        cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 {grid_step} {grid_n} 104\n"
        #rc = blocking_command(cfg,cmd)
        system_command(cfg, cmd)


def mp_rf1(cfg,multis, ras, decs,num_procs=NumProcs_mp_rf1):
    """

    these are running basically blind ... the output is to files and I only need to know that it is done,
    regardless of success or fail

    so, break up in to XXX tasks here and spin up, then return when all done

    :param cfg:
    :return:
    """
    #cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 104\n"
    #system_command(cfg, cmd)

    #num_jobs = np.count_nonzero(len(multis))
    #num_multi_per_job = int(np.ceil(num_jobs / num_threads))  #min(max_threads, num_jobs // nominal_jobs_per_process)

    #print(f"Top of mp_rf1: len(multis) = {len(multis)}")

    #tune num procs
    active_shots = node_active_ct(cfg)
    if active_shots > 0:
        num_procs = int(np.floor(MaxTotalProcs_mp_rf1 / active_shots))
        num_procs = max(1,min(num_procs, MaxPerShotProcs_mp_rf1)) #always at least one

    print(f"[{cfg.datevshot}] line detections ***MULTITHREADED*** (rdet_rf1) ({len(ras)}), using num_procs = {num_procs}...")

    with multiprocessing.Manager() as manager:
        mgr_list = manager.list()
        idx = np.array_split(np.arange(len(multis)),num_procs)
        processes = []

        for i in range(num_procs):
            print(f"{cfg.datevshot} (mp_rf1) : Spinning up rf1 i={i},  len(multis) = {len(multis[idx[i]])}, idx[i] = {idx[i]}")
            process = multiprocessing.Process(target=mp_rf1_worker,
                                               args=(mgr_list,cfg,i,idx[i],multis[idx[i]],ras[idx[i]],decs[idx[i]]))

            processes.append(process)
            process.start()

        # Wait for processes to complete
        for process in processes:
            print(f"[{cfg.datevshot}] (mp_rf1) : Joining {process}")
            process.join()

        print(f"[{cfg.datevshot}] (mp_rf1) : rdet_rf1 all done")
        #again, don't care about the results here, just need them done

def rdet_rf1(cfg):
    """
    emission line detections

    :param cfg:
    :return:
    """

    try:
        #rf1 244.7716830 33.8165741 35 4505 50 513_105_051 20240730v009 1.70 3.0 3.5 0.5 3 104
        #rf1 244.7050930 33.8114395 35 4505 50 514_103_019 20240730v009 1.70 3.0 3.5 0.5 3 104

        failed_list = []
        #passed_list = []
        output_extensions = [".list", ".mc", ".spec"]

        rc = 0
        os.chdir(os.path.join(cfg.cwd, "alldet"))  # make sure we are in the right directory
        os.makedirs("detect_out", exist_ok=True)

        multis = np.loadtxt(os.path.join("../getcen/",f"ifucen_{cfg.datevshot}.dat"),dtype=str,usecols=(0),unpack=True)
        ras,decs = np.loadtxt(os.path.join("../getcen/", f"ifucen_{cfg.datevshot}.dat"), dtype=float,usecols=(1,2),unpack=True)

        # print(f"line detections (rdet_rf1) ({len(ras)})...")
        # #todo: turn this next loop into a multithread call
        # ct = 0
        # for multi, ra, dec in zip(multis, ras, decs):
        #     ct +=1
        #     print(f"{ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ")#,end="")
        #     cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 104\n"
        #     system_command(cfg,cmd)
        #
        #     # check the output exists cal_out/20240730v009_514_103_019_cal.fits
        #
        #     #todo: assume good for now ... need to check
        #     #print("done (todo: check output files)")
        #     #check that .list, .mc, .spec exist
        #     #20240730v009_025_067_032.list

        #print(f"[{cfg.datevshot}] line detections ***MULTITHREADED*** (rdet_rf1) ({len(ras)})")#, num_procs = {NumProcs_mp_rf1}...")
        mp_rf1(cfg,multis, ras, decs)#,num_procs=6)

        #grid_n = 3  # n x n so 3 = a 3x3 grid
        #grid_step = 0.5  # grid step size
        #print(f"[{cfg.datevshot}] *** test *** set grid_step and grid_n to 0.0 1")
        #DD 20260627 ... we like 1 and 0, so keep that as the new default
        # grid_n = 1 , grid_step = 0.0
        grid_n = f"{int(cfg.linedet_parms[0])}"
        grid_step = f"{float(cfg.linedet_parms[1]):0.1f}"

        # if cfg.numexp < 3 or cfg.dither_configuration == 0:
        #     grid_n = 13
        #     grid_step = 0.25

        #now - check them all (this can be serial, it is fast)
        ct = 0
        for multi, ra, dec in zip(multis, ras, decs):
            ct +=1
            base_str = f"[{cfg.datevshot}] checking rdet_rf1 ({ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... "
            #print(f"{ct}) checking [{cfg.datevshot}] {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ", end="")
            output_found = np.array([0,0,0])
            for i,ext in enumerate(output_extensions):
                outfn = f"{cfg.datevshot}_{multi[6:]}{ext}"

                if os.path.exists(os.path.join("detect_out/", outfn)):
                    output_found[i] = 1

            if np.count_nonzero(output_found) != 3:
                #something failed, we will want to re-run these once
                cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 {grid_step} {grid_n} 104\n"
                failed_list.append(cmd)
                print(f"{base_str} FAIL. May re-run at the end.")
            else:
                print(f"{base_str} pass")

        #let any that failed re-run as serial to keep it simple
        if len(failed_list) > 0:
            if len(failed_list) < len(ras):
                #print(f"{base_str} {len(failed_list)} failed. Can be transient issues, so will re-run ...")
                print(f"{cfg.datevshot} (rdet_rf1) : serially re-run {len(failed_list)} failed IFUs ...")
                #reset rc
                rc = 0
                for ct,cmd in enumerate(failed_list):
                    print(f"{cfg.datevshot} (rdet_rf1) : {ct} {cmd.split()[1]} {cmd.split()[2]} {cmd.split()[6]} ...")#,end="")
                    system_command(cfg, cmd)

                    output_found = np.array([0, 0, 0])
                    for i, ext in enumerate(output_extensions):
                        outfn = f"{cfg.datevshot}_{cmd.split()[6]}{ext}"
                        if os.path.exists(os.path.join("detect_out/", outfn)):
                            output_found[i] = 1

                    if np.count_nonzero(output_found) != 3:
                        # something failed, we will want to re-run these once
                        #failed_list.append(cmd)
                        print(f"{cfg.datevshot} (rdet_rf1) : {ct} {cmd.split()[1]} {cmd.split()[2]} {cmd.split()[6]} ... FAIL. Second attempt. No more retries.")
                        rc = 1 #some failures, but not all
                    else:
                        print(f" pass")
            else:
                print(f"[{cfg.datevshot}] (rdet_rf1) All failed. Will not attempt full re-run.")

    except:
        print(traceback.format_exc())
        rc = -1
        print(f"{cfg.datevshot} (rf1) : Fatal exception in rdet_rf1!")

    return rc

def approx_snr(snr_array):
        """
        a curve that approximates a nominal ***MAXIMUM*** of what we expect from a normal SNR
            (normally from 4.8 to 10.0) for a single IFU, assuming 3x 360sec exposures and good seeing, etc.
            note that there can be factors of 2x-3x just by varying the conditions a bit
            and, of course, any REAL line of sight variations

        imposed conditions:

        (continuum around line > -3)   &
        (4.8 <= snr <= 10.0) &
        (chi2 < 1.5 ) &
        (1.5 <= linewidth(sigma) <= 16.) &
        (chi2_fib <= 4.5)

        ignored conditions: seeing, tput, etc (just assume to be good)

        :param snr_array: (snr array, normally 4.8 to 10.0 in steps of 0.1)
        :return: counts expected per IFU (under normal 3-dither, 360second exposures)
        """

        #the np.exp(-0.9 * SNR) shape is pretty good (at SNR << 4.5, the upturn should get sharper though)
        #the scale in front could be 150 to 450, depending on number of factors (See above note)
        #so, run with 300. , in the middle, for now
        return 300. * np.exp(-0.9 * snr_array) #per IFUs though, down to snr 4.8 with limited validation


def approx_count_at_snr(snr_array, total_exp_time, num_dithers):
        """

        #adjust the SNR based on approximate proportionality of sqrt(time)
        #e.g. adjust the depth in a way, based off the empirical model and the standard HETDEX 360 second exposure
        # since we are using the TOTAL time here, need to set that as per exposure
        # AND we assume, if multiple exposuress that they are dithered and do not overlap

        :param snr_array:
        :param total_exp_time:
        :param num_dithers:
        :return: adjusted SNRs and counts expected per IFU
        """

        per_exp_time_mux = np.sqrt(total_exp_time / num_dithers / 360.0)

        # adjusts the area covered
        vol_mux = num_dithers / 3

        #shift the SNR array by 1/srt(time) and multiply by volume change
        return approx_snr(snr_array/per_exp_time_mux) * vol_mux


#
# defunct? replaced with  check_line_detections_by_ifu
#
# def check_line_detections(cfg):
#     """
#     basic examination of results of line detecetions (rdet_rf1)
#
#     summary --- if there are too many line detections (with adjustments for exposure time, SNR, etc)
#                 then either warn and continue or fail and abort
#
#     not worried about too few detections at this time
#
#     :param cfg:
#     :return:
#     """
#
#     #todo: switch to per IFU
#     #      use biweight instead of sum ... so we need to track the dets by IFU to do this
#
#
#
#     rc = 0
#     min_snr = 4.8
#     max_snr = 10.0
#     step_snr = 0.1
#     warn_thresh = 3.0
#     fail_thresh = 10.0
#     try:
#         #recall, .mc has the 1 line per detction (.spec has the spectra, .list has all the involved fibers)
#         #        colnames = ['wave', 'wave_err','flux','flux_err','linewidth','linewidth_err',
#         #             'continuum','continuum_err','sn','sn_err','chi2','chi2_err','ra','dec',
#         #             'datevshot','noise_ratio','linewidth_fix','chi2_fix', 'chi2fib',
#         #             'src_index','multiname', 'exp','xifu','yifu','xraw','yraw','weight',
#         #             'apcor','sn_cen', 'flux_noise_1sigma', 'sn_3fib', 'sn_3fib_cen','dummy']
#
#         fns = glob.glob(os.path.join(cfg.cwd, "alldet/detect_out/*.mc"))
#
#         #one per IFU
#         num_ifus = len(fns)
#         all_lw = []
#         all_cont = []
#         all_snr = []
#         all_chi2 = []
#         all_chi2_fib = []
#
#         for fn in fns:
#             xlw, xcont, xsnr, xchi2, xchi2fib = np.loadtxt(fn, usecols=[4, 6, 8, 10, 18], unpack=True)
#             all_lw += list(xlw)
#             all_cont += list(xcont)
#             all_snr += list(xsnr)
#             all_chi2 += list(xchi2)
#             all_chi2_fib += list(xchi2fib)
#
#         all_lw = np.array(all_lw)
#         all_cont = np.array(all_cont)
#         all_snr = np.array(all_snr)
#         all_chi2 = np.array(all_chi2)
#         all_chi2_fib = np.array(all_chi2_fib)
#
#         # default: based on snr >=4.8
#         #all_snr = np.array(sorted(all_snr[all_snr >= min_snr]))
#
#         #basic, standard selection, but not exhaustive
#         #!!! make sure any changed match with check_line_detections_by_ifu() !!!
#         sel = np.array(all_cont > -3.) * np.array(all_snr >= 4.8) * np.array(all_snr <= 10.0) * np.array(all_chi2 < 1.5) * \
#               np.array(all_lw >= 1.5) * np.array(all_lw <= 16.) * np.array(all_chi2_fib <= 4.5)
#         all_snr = np.array(sorted(all_snr[sel]))
#
#         if len(all_snr) > 0:
#             data_min_snr = round(np.nanmin(all_snr),1)
#             data_max_snr = round(np.nanmax(all_snr),1)
#
#             min_snr = max(min_snr,data_min_snr)
#             max_snr = min(max_snr, data_max_snr)
#             xbins = np.arange(min_snr,max_snr+step_snr,step_snr)
#
#             binned_snr = np.histogram(all_snr,bins=xbins)
#             data_snr = np.sum(binned_snr[0])
#
#             #model_snr = int(np.sum(approx_snr(xbins)) * np.sqrt(cfg.total_exp_time / 1080.) * num_ifus / 78)
#
#             model_snr = int(np.sum(approx_count_at_snr(xbins,cfg.total_exp_time,cfg.numexp))) * num_ifus
#             print(f"[{cfg.datevshot}] Line dets for SNR inputs: IFUs = {num_ifus}, exp time = {cfg.total_exp_time},"
#                   f"num exp = {cfg.numexp}, SNR range = [{xbins[0]},{xbins[-1]}], result = {model_snr}")
#             print(f"[{cfg.datevshot}] Line dets for SNR [{min_snr:0.1f},{max_snr:0.1f}]. Data / Model {data_snr} / {model_snr} ")
#
#             cfg.ratio_line_dets = data_snr / model_snr
#
#             if cfg.ratio_line_dets > fail_thresh:
#                 if not FORCE_CONTINUE:
#                     print(f"[{cfg.datevshot}] Fail! {cfg.ratio_line_dets:0.1f}x number of expected detections "
#                           f"SNR [{min_snr:0.1f},{max_snr:0.1f}]. Will terminate.")
#                     rc = -1
#                 else:
#                     print(f"[{cfg.datevshot}] Warning! {cfg.ratio_line_dets:0.1f}x number of expected detections "
#                           f"SNR [{min_snr:0.1f},{max_snr:0.1f}]. "
#                           f"Force flag set, so will continue")
#                     rc = 1
#             elif cfg.ratio_line_dets > warn_thresh:
#                 print(f"[{cfg.datevshot}] Warning! {cfg.ratio_line_dets:0.1f}x number of expected detections "
#                       f"SNR [{min_snr:0.1f},{max_snr:0.1f}]")
#                 rc = 1
#             else: # probably fine here
#                 print(f"[{cfg.datevshot}] (DEBUG) {cfg.ratio_line_dets:0.1f}x number of expected detections "
#                       f"SNR [{min_snr:0.1f},{max_snr:0.1f}]")
#                 rc = 0
#
#             print(f"[{cfg.datevshot}] !!! DEBUG !!! for now, force the rc to be zero")
#             rc  =0
#
#     except:
#         print(traceback.format_exc())
#         rc = 0 #cannot make a call here
#
#     return rc


def adjust_minimum_snr(single_exposure_time, avg_sky, seeing = None,
                       baseline_snr=4.8, baseline_sky=500., baseline_seeing= 1.8):
    """

    based on the (approximate) depth as scaled by sqrt(single exposure time / 360.0s baseline)
    and under the assumption that this holds and shot (sky) noise does not rise faster ...
       we do, however, attempt to correct for extra sky noise with a term that divides by the
       sqrt of the ratio of the shot average sky / baseline sky

    the idea is that for a normal HETDEX exposure of about 360s (1080 sec total for 3 dithers),
      SNR ~ 4.8 is about the lowest we care to go

    for longer exposures, assuming the sqrt(time) holds, those that would have been 4.8 (real) detections
      get larger SNR. For 2250 single expsoure (6.25x longer), that is like 4.8 * sqrt(2250/360s) ~ 12.0 SNR
      (note there are exposures up into the 4000-5000s range)
      BUT this has, say, avg_sky around 1750, so 12.0 SNR * 1/sqrt(1750/400)) ~ 5.7 SNR

    and stuff we would consider junk/noise can rise to the minimum very easily

    --- OR --- INVERT this and get a minimum SNR that knocks out junk??

    HMMM Maybe we need to INCREASE rather than decress when sky is bad? the idea being that bad sky
     makes it easier to get a higher SNR that is just noise? so we'd raise the SNR floor to cut those out?

     so longer exposures means deeper and maybe we keep SNR the same, BUT bad sky increases noise more
       so we raise the floor?


    DOES not consider seeing, tput or apcor

    :param single_exposure_time: the exposure time for a SINGLE exposure (dither)
    :param avg_sky: the average sky for the shot
    :param baseline_snr: the SNR (e.g. the normal floor)
    :param baseline_sky: a baseline reference sky value (most HETDEX are in the low 100s)
    :return:
    """

    #these two conteract each other ...
    # with time_adjust allowing for what would be a lower SNR for a shorter exposure (signal increase) and
    # sky_adjust pushing that SNR higher when the sky is worse (noise increase)
    # i.e. for a long exposure, if the sky remained low, you could start with lower SNR which would be the equivalent
    #   of 4.8 on a standard exposure, but for large enough sky, regardless of exposure length, the sky noise is
    #   high enough to generate too many low SNR (false) detections


    # todo --OR-- assume signal naturally increases and we don't need to adjust by time
    #             AND we only correct UP (make the lower cutoff a larger SNR) based on the avg_sky
    #                       with a higher baseline_sky (maybe around 1000)
    #    note: HETDEX length exposures have avg_sky 120 - 330 range
    # if we dump time_adjust and just use sky, then a baseline around 1000 is about right?

    sky_adjust = np.sqrt(avg_sky/baseline_sky)
    time_adjust = np.sqrt(single_exposure_time / 360.0) #for long expsousres adjust target down
    if seeing is not None: #for good seeing this is < 1 and adjusts the target SNR down
        seeing_adjust = seeing / baseline_seeing
    else:
        seeing_adjust = 1.0
    net_adjust = sky_adjust / time_adjust

    return max(baseline_snr, np.round(baseline_snr * net_adjust * seeing_adjust,1))
    #return np.round(baseline_snr * net_adjust, 1)

def snr_rescale(avg_sky,baseline_sky=1000.):
    """
    simple scaling based on the avg_sky

    this would be used to re-scale the *.mc file reported SNR
    for really bad sky, the idea is that the SNR get much smaller
    (e.g. for sky around 7500.0 and a baseline of 1000.0, reported SNR 13 becomed about 4.8)
    empirically, so far, this seems about right (with say 50%?) that is it might be for sky 7500.0 we
      want something between SNR 13 and 17 to go to 4.8

    so you'd take this and multiply it through the *.mc files ... need to feed that to HETDEX_API (otherwise
      it takes way to long to build up the line.h5 file)
        #>>> since we can pass in an SNR cut, might flip this around and pass in the SNR that WOULD be 4.8
        #>>> and then rescale the resulting values? ... NOTE: if ELiXer runs these, it would get the higher
        #>>> SNR, I think ... could instead raise the flux_err  arrays?

    the seeing and the tput and the exposure length take care of themselves, boosting the signal to noise

    this is a correction to attempt to handle un-accounted for really bad sky "noise" (not a sky residual, per se)

    :param avg_sky:
    :param baseline_sky:
    :return:
    """

    #as showm, this returns a value >= 1.0
    #with the intended use being to divide the *.mc SNR columns to lower their value and keep a fixed
    #  snr minimum or multiply into the snr minimum to raise the lower cut off
    return max(1.0, np.sqrt(avg_sky / baseline_sky))

#def raise_snr_floor(single_exposure_time, avg_sky, baseline_snr=4.8, baseline_sky=1000.):

def get_raw_line_detections_from_mc(cfg):
    """

    notice! a version of this is also in check_line_detctions_by_ifu() since it
            also iterates over the mc files

    :param cfg:
    :return: total # of line detects, line detects by bin, left edges of the bins
    """

    step = 0.1
    counts = []
    edges = []
    total = -1
    try:
        mc_path = os.path.join(cfg.cwd_orig,f"sci{cfg.datevshot}/alldet/detect_out")
        if os.path.exists(mc_path):
            fns = glob.glob(f"{mc_path}/*.mc")
            all_sn = []
            for fn in fns:
                sn = np.loadtxt(fn, usecols=[8])
                all_sn += list(sn)
            all_sn = np.array(all_sn)

            min_sn = np.min(all_sn)
            max_sn = np.max(all_sn)
            bin_sn = np.arange(np.floor(min_sn * 10.)/10.,np.ceil(max_sn * 10.)/10.+step,step)
            counts, edges = np.histogram(all_sn,bins=bin_sn)

            if len(edges) > 0:
                edges = edges[:-1] #cut off the right most edge so counts matches (left) edges
    except:
        print(traceback.format_exc())

    return total, counts, edges

def check_line_detections_by_ifu(cfg):
    """
    basic examination of results of line detecetions (rdet_rf1) by IFU
    (same as check_line_detection) but rather than as the whole shot, is IFU by IFU

    summary --- if there are too many line detections (with adjustments for exposure time, SNR, etc)
                then either warn and continue or fail and abort

    not worried about too few detections at this time

    :param cfg:
    :return:
    """


    rc = 0
    min_snr = 4.8
    max_snr = 10.0
    step_snr = 0.1
    warn_thresh = 5.0
    fail_thresh = 10.0
    try:
        #recall, .mc has the 1 line per detction (.spec has the spectra, .list has all the involved fibers)
        #        colnames = ['wave', 'wave_err','flux','flux_err','linewidth','linewidth_err',
        #             'continuum','continuum_err','sn','sn_err','chi2','chi2_err','ra','dec',
        #             'datevshot','noise_ratio','linewidth_fix','chi2_fix', 'chi2fib',
        #             'src_index','multiname', 'exp','xifu','yifu','xraw','yraw','weight',
        #             'apcor','sn_cen', 'flux_noise_1sigma', 'sn_3fib', 'sn_3fib_cen','dummy']

        fns = glob.glob(os.path.join(cfg.cwd, "alldet/detect_out/*.mc"))

        #one per IFU
        num_ifus = len(fns)

        data_ct = [] #one per IFU
        model_ct = []

        all_sn = []



        for fn in fns:
            xlw, xcont, xsnr, xchi2, xchi2fib = np.loadtxt(fn, usecols=[4, 6, 8, 10, 18], unpack=True)
            ifu_name = os.path.basename(fn).split(".")[0].split("_")[1:] #ie. 20240801v002_323_043_040.mc
            ifu_name = '_'.join(ifu_name)
            cfg.ifu_list.append(ifu_name)
            cfg.ifu_linedet_ct.append(-1)
            cfg.ifu_linedet_ratio.append(-1)

            all_sn += list(xsnr)


            #this is a basic sanity sub-selection ... it may be different than what is applied
            #for the input to ELiXer, but, this is a way to semi-calibrate these data
            #to what is rouhgly, normally expected to see if we are way off
            # (i.e. we would normally expecte something like 5-20 line detects per IFU to survive this
            #       down-selection for a 3-Dither, 1100s observation)

            ## !!!! make sure any changed MATCH those in check_line_detections() !!!
            sel = np.array(xcont > -3.) * np.array(xsnr >= 4.8) * np.array(xsnr <= 10.0) * np.array(xchi2 < 1.5) * \
                  np.array(xlw >= 1.5) * np.array(xlw <= 16.) * np.array(xchi2fib <= 4.5)

            xsnr = sorted(xsnr[sel])

            if len(xsnr) > 1:
                data_min_snr = round(np.nanmin(xsnr),1)
                data_max_snr = round(np.nanmax(xsnr),1)

                bin_min_snr = max(min_snr,data_min_snr)
                bin_max_snr = min(max_snr, data_max_snr)
                xbins = np.arange(bin_min_snr,bin_max_snr+step_snr,step_snr)

                #does not really have to be a histogram, but we might want to see this in the future
                #plus this handles summing just over the desired SNR range
                binned_snr = np.histogram(xsnr,bins=xbins)
                data_ct.append(np.sum(binned_snr[0]))
                cfg.ifu_linedet_ct[-1] = data_ct[-1]

                # build a model, based on the SNR range in the IFU (each can be different)
                # model_snr = int(np.sum(approx_snr(xbins)) * np.sqrt(cfg.total_exp_time / 1080.))
                model_snr = max(1,int(np.sum(approx_count_at_snr(xbins, cfg.total_exp_time, cfg.numexp))))
                cfg.ifu_linedet_ratio[-1] = data_ct[-1] / max(1,model_snr)

            else: #NO appropriate data in IFU
                #zero_ifu.append(os.path.basename(fn))
                print(f"[{cfg.datevshot}] ***DEBUG*** No detections in range for IFU: {os.path.basename(fn)}")
                data_ct.append(len(xsnr)) #either 0 or 1 here
                cfg.ifu_linedet_ct[-1] = data_ct[-1]
                #xbins = np.arange(min_snr, max_snr + step_snr, step_snr)
                model_snr = 1 #just so we won't get a divide by zero, but the real count will be 0 or 1 anyway
                cfg.ifu_linedet_ratio[-1] = 0

            model_ct.append(model_snr)

        all_sn = np.array(all_sn)

        min_mc_snr = np.min(all_sn)
        max_mc_snr = np.max(all_sn)
        bin_sn = np.arange(np.floor(min_mc_snr * 10.) / 10., np.ceil(max_mc_snr * 10.) / 10. + 0.1, 0.1)
        counts, edges = np.histogram(all_sn, bins=bin_sn)
        total_mc_line_dets = len(all_sn)

        #this will be a stupid long string
        snr_log_str = ", ".join([f"{e:0.1f} ({c})" for e, c in zip(edges, counts)])
        print(f"[{cfg.datevshot}] Raw *.mc line detects by snr. Total = {total_mc_line_dets} for "
              f"SNR {min_mc_snr:0.1f} to {max_mc_snr:0.1f}: {snr_log_str}")

        #todo: make a suggested SNR cut based on the avg_sky?? or other parameters ...
        #      research on that formulation still in progress
        if cfg.avg_sky is None:
            cfg.avg_sky = get_avg_sky(cfg)
        cfg.snr_rescale = snr_rescale(cfg.avg_sky)

        if cfg.snr_rescale > 1.0:
            print(f"[{cfg.datevshot}] Recommend rescale minimum SNR. "
                  f"Line count @ 4.8 = {np.count_nonzero(all_sn>=4.8)}. "
                  f"Line count @ {4.8 * cfg.snr_rescale:0.1f} = {np.count_nonzero(all_sn>=(4.8*cfg.snr_rescale))}")

        if len(edges) > 0:
            edges = edges[:-1]  # cut off the right most edge so counts matches (left) edges

        #since the SNR range varies from IFU to IFU, cannot report a single range
        #NOR is it correct to just run one estimate and multiply by the number of IFUs
        print(f"[{cfg.datevshot}] Model line dets for SNR inputs (nominal {min_snr:0.1f} to {max_snr:0.1f}): "
              f"IFUs = {num_ifus}, exp time = {cfg.total_exp_time:0.1f},"
              f"num exp = {cfg.numexp}, total count = {np.sum(model_ct)}")

        #now compare data vs model
        data_ct = np.array(data_ct)
        model_ct = np.array(model_ct)

        sel = model_ct > 0 #should be all of them
        data_ct = data_ct[sel]
        model_ct = model_ct[sel]

        data_over_model = data_ct / model_ct
        cfg.ratio_line_dets = np.sum(data_ct) / np.sum(model_ct) #total ratio

        fail_ct = np.count_nonzero(data_over_model >= fail_thresh)
        warn_ct = np.count_nonzero(data_over_model >= warn_thresh) - fail_ct
        pass_ct = np.count_nonzero(data_over_model < warn_thresh)

        print(f"[{cfg.datevshot}] Downselected line det counts: totals data = {np.sum(data_ct)} vs model = {np.sum(model_ct)}; "
              f" by IFU: {fail_ct} Fail, {warn_ct} Warn, {pass_ct} Pass")

        if fail_ct > 0 or warn_ct > 0:
            idx_ifu = np.where(data_over_model >= warn_thresh)[0]
            for i in idx_ifu:
                stat_str = "fail" if cfg.ifu_linedet_ratio[i] >= fail_thresh else "warn"
                print(f"[{cfg.datevshot}] Line Dets {stat_str}: IFU {cfg.ifu_list[i]} has {cfg.ifu_linedet_ct[i]} dets at "
                      f"{cfg.ifu_linedet_ratio[i]}x the maximum expected.")

        #todo: could we flag an IFU as bad for this reason? That its line dets are way too high?
        # note: amp_stats() have already been built ... could append to it?

        #if most are fails or warns, fail it. Let the fails count twice when combining with warns
        if (fail_ct / num_ifus >= 0.5) or \
                ( (fail_ct / num_ifus >= 0.2) and (2 * fail_ct + warn_ct) / num_ifus >= 1.0):
            if not FORCE_CONTINUE:
                print(f"[{cfg.datevshot}] Fail! {cfg.ratio_line_dets:0.1f}x nominal maximum total number of"
                      f" expected detections SNR [{min_snr:0.1f},{max_snr:0.1f}]. Will terminate.")
                rc = -1
            else:
                print(f"[{cfg.datevshot}] Warning! {cfg.ratio_line_dets:0.1f}x nominal maximum total number of"
                      f" expected detections SNR [{min_snr:0.1f},{max_snr:0.1f}]. "
                      f"--force flag set, so will continue")
                rc = 1
        elif warn_ct / num_ifus >= 0.5: #if half or more are warns
            print(f"[{cfg.datevshot}] Warning! {cfg.ratio_line_dets:0.1f}x nominal maximum total number of"
                  f" expected detections SNR [{min_snr:0.1f},{max_snr:0.1f}]")
            rc = 1
        else:
            print(f"[{cfg.datevshot}] (DEBUG) {cfg.ratio_line_dets:0.1f}x nominal maximum total number of"
                  f" expected detections")
            rc = 0

    except:
        print(traceback.format_exc())
        rc = 0 #cannot make a call here

    return rc



#continuum detection
def rgetmax(cfg):
    """
    continuum detections

    the datevsshot/cs/rgetmax script changes directories and runs under /tmp/maxflux<datevshot>
      each exp re-uses the same directory
      then all exps are rolled up

    :param cfg:
    :return:
    """

    try:
        #failed_list = []
        #passed_list = []
        output_extensions = [".list", ".mc", ".spec"]

        rc = 0
        os.chdir(os.path.join(cfg.cwd, "cs"))  # make sure we are in the right directory
        os.makedirs("spec", exist_ok=True)

        print(f"[{cfg.datevshot}] continuum detections (rgetmax) ...")

        mux = linear_exptime_scale(cfg)
        if mux > 1.0:
            mux = np.sqrt(mux)
            print(f"[{cfg.datevshot}] updating for extended exposure time by x{mux:0.1f}")
            cutstr = f"cut={20.0 * mux:0.1f}"
            system_command(cfg, f"sed -i s/cut=20./{cutstr}/ rgetmax")


        if cfg.exp == 0 and cfg.numexp == 3:
            print(f"[{cfg.datevshot}] using standard 3-dither rgetmax")
            #pass #do nothing, the default of 3 exposures stands
        else:
            #note: technically, this could fail if this is on a --resume and the original
            # string was already replaced, but in that case, the sed just won't do anything
            # and we would already have the correct exposures to (re)run
            print(f"[{cfg.datevshot}] updating for selected and available exposures")
            #exp_array = ("exp01" "exp02" "exp03")
            #sed -i s#exp_array=\(\"exp01\"\ \"exp02\"\ \"exp03\"\)#exp_array=\(\"exp01\"\)# rgetmax
            cmd = f"sed -i s#exp_array=\(\\\"exp01\\\"\ \\\"exp02\\\"\ \\\"exp03\\\"\)#"
            if cfg.exp > 0:
                cmd += f"exp_array=\(\\\"exp{str(cfg.exp).zfill(2)}\\\"\)# rgetmax"
            else:
                cmd += f"exp_array=\("
                for exp in range(cfg.numexp):
                    cmd += f"\\\"exp{str(exp+1).zfill(2)}\\\"\ "
                cmd += f"\)# rgetmax"

            system_command(cfg,cmd)

        print(f"[{cfg.datevshot}] running rgetmax (continuum detection) ...")

        system_command(cfg,f"rgetmax {cfg.datevshot.split('v')[0]} {cfg.datevshot.split('v')[1]}")

        #output is a cs.tar file ... not broken up by IFU #20240730v009cs.tar
        if os.path.exists(f"spec/{cfg.datevshot}cs.tar"):
            #assume good
            #contains 1 or more .spec and .list files (one per detection)
            #a single .cs file (may be coords and counts ... but is longer than expected)
            #as well as a single .rcs file that seems to just contain excecutable calls to rf1 (prob. to line extract at the position)

            #need to untar for use
            system_command(cfg, f"cd spec; tar -xf {cfg.datevshot}cs.tar ; cd ..")
        else:
            #a failure, but not fatal
            rc = 1
            print(f"[{cfg.datevshot}] Continuum detections fail. Not fatal. cs/spec/{cfg.datevshot}cs.tar not created.")


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Fatal exception in rgetmax.",traceback.format_exc())

    return rc



def build_shot_h5(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        print(f"[{cfg.datevshot}] Constructing shot hdf5 file. This will take a while ... ")
        rc = 0

        #some setup that hetdex_api needs
        os.chdir(os.path.join(cfg.cwd,"vdrp"))

        cmd = f"ln -s shifts/{cfg.starcat_ast}/{cfg.datevshot}/shout.* ."
        system_command(cfg, cmd)

        cmd = f"ln -s shifts/{cfg.starcat_ast}/{cfg.datevshot}/2* ."
        system_command(cfg, cmd)

        cmd = f"ln -s shifts/{cfg.starcat_ast}/{cfg.datevshot}/all.mch ."
        system_command(cfg, cmd)

        cmd = f"ln -s shifts/{cfg.starcat_ast}/{cfg.datevshot}/radec* ."
        system_command(cfg, cmd)

        cmd = f"ln -s ../detect/{cfg.datevshot}/norm.dat ."
        system_command(cfg, cmd)

        os.makedirs("match_pngs", exist_ok=True)
        os.chdir(os.path.join(cfg.cwd, "vdrp/match_pngs"))
        cmd = f"ln -s ../shifts/{cfg.starcat_ast}/{cfg.datevshot}/match* ."
        system_command(cfg, cmd)

        os.chdir(cfg.cwd) #also need a match_pngs up one level
        os.makedirs("match_pngs", exist_ok=True)
        os.chdir(os.path.join(cfg.cwd, "match_pngs"))
        cmd = f"ln -s ../vdrp/shifts/{cfg.starcat_ast}/{cfg.datevshot}/match* ."
        system_command(cfg, cmd)

        # print("Copying match_pngs .... ")
        # try:
        #     #note: just copy the pdfs ... hetdex_api create_astrometry handles the conversion to png and renaming
        #     done = False
        #     #GAIA first
        #     mp_path = os.path.join(cfg.cwd,f"vdrp/shifts/gaia/{cfg.datevshot}")
        #     if os.path.exists(mp_path):
        #         match_pdfs =  glob.glob(f"{mp_path}/match_exp*.pdf")
        #         if len(match_pdfs) > 0:
        #             system_command(cfg, f"cp {mp_path}/match_exp*.pdf {cfg.cwd}/match_pngs")
        #             done = True
        #             print("Used GAIA match_pngs")
        #
        #     if not done:
        #         # SDSS next
        #         mp_path = os.path.join(cfg.cwd, f"vdrp/shifts/sdss/{cfg.datevshot}")
        #         if os.path.exists(mp_path):
        #             match_pdfs = glob.glob(f"{mp_path}/match_exp*.pdf")
        #             if len(match_pdfs) > 0:
        #                 system_command(cfg, f"cp {mp_path}/match_exp*.pdf {cfg.cwd}/match_pngs")
        #                 done = True
        #                 print("Used SDSS match_pngs")
        #
        #
        #     if not done:# PanSTARRS last
        #         mp_path = os.path.join(cfg.cwd, f"vdrp/shifts/panstarrs/{cfg.datevshot}")
        #         if os.path.exists(mp_path):
        #             match_pdfs = glob.glob(f"{mp_path}/match_exp*.pdf")
        #             if len(match_pdfs) > 0:
        #                 system_command(cfg, f"cp {mp_path}/match_exp*.pdf {cfg.cwd}/match_pngs")
        #                 done = True
        #                 print("Used PanSTARRS match_pngs")
        #
        #     if not done:
        #         print("WARNING!!! Non-fatal. Could not find match_pngs under vdrp/shifts/<catalog>/<datevshot>")
        # except:
        #     print(traceback.format_exc())


        #########################
        # initial hdf5 file
        ########################

        os.chdir(cfg.cwd)

        cmd = f"python3 {hetdex_api_path}/h5tools/create_shot_hdf5.py"
        cmd += " --tar"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" --rootdir \"{cfg.cwd}\""

        if cfg.starcat_ast =='gaia' and os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.gaia"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.gaia\""
            cfg.starcat_ast = "gaia"
        elif cfg.starcat_ast =='sdss' and os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.sdss"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.sdss\""
            cfg.starcat_ast = "sdss"
        elif cfg.starcat_ast =='panstarrs' and os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.panstarrs"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.panstarrs\""
            cfg.starcat_ast = "panstarrs"
        else: #should not happen, but choose one, if exists
            if  os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.gaia"):
                cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.gaia\""
                cfg.starcat_ast = "gaia"
            elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.sdss"):
                cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.sdss\""
                cfg.starcat_ast = "sdss"
            elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.panstarrs"):
                cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.panstarrs\""
                cfg.starcat_ast = "panstarrs"
            else:
                print("Fatal: cannot find *.dithall file for this shot")
                return -1

        system_command(cfg, cmd)

        #assume good?
        if not os.path.exists(f"{cfg.datevshot}.h5"):
            rc = -1
            return rc

        print(f"[{cfg.datevshot}] Created: {cfg.cwd}/{cfg.datevshot}.h5")

        #########################
        # now append_calfib
        ########################
        print(f"[{cfg.datevshot}] Appending calibrated fibers ... ")
        cmd = f"python3 {hetdex_api_path}/h5tools/append_calfib.py"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" --rootdir \"{cfg.cwd}/alldet/cal_out/\""

        system_command(cfg, cmd)

        #assume okay?


        #########################
        # now create_fullsky_model
        ########################
        print(f"[{cfg.datevshot}] Appending fullsky model ... ")
        cmd = f"python3 {hetdex_api_path}/h5tools/create_fullskymodel_hdf5.py"
        cmd += " --append"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" -r \"{cfg.cwd}/alldet/output/\""

        system_command(cfg, cmd)


        #########################
        # now create_cal_hdf5
        ########################
        print(f"[{cfg.datevshot}] Create cal hdf5  ... ")


        #hetdex_api needs the local fwhm.all in a different format
        #detect/20240730v009
        if not os.path.exists(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail"):
            system_command(cfg,f"mv {cfg.cwd}/detect/{cfg.datevshot}/fwhm.all {cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail")
        #only reads the first row
        fwhm, err, ns = np.loadtxt(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail",unpack=True,usecols=[0,1,2],
                                   max_rows=1,dtype=float)
        try:
            with open(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.all","w") as f:
                f.write(f"{cfg.datevshot} {fwhm} {err} {int(ns)}\n")
        except:
            print(f"[{cfg.datevshot}] Fatal. Could not build fwhm.all file from fwhm.detail. fhwm, err, ns = {fwhm} {err} {ns}.")
            #print(traceback.format_exc())
            rc = -1
            print(f"[{cfg.datevshot}] Could not build hdf5 shot file.",traceback.format_exc())
            return rc


        cmd = f"python3 {hetdex_api_path}/h5tools/create_cal_hdf5.py"
        cmd += " --append"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" -tp \"{cfg.cwd}/detect/\""
        cmd += f" -detdir  \"{cfg.cwd}/detect/{cfg.datevshot}/\""

        system_command(cfg, cmd)



        #########################
        # now create_astrometry_hdf5
        ########################
        print(f"[{cfg.datevshot}] Create astrometry  ... ")

        cmd = f"python3 {hetdex_api_path}/h5tools/create_astrometry_hdf5.py"
        cmd += " --append"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" -detdir  \"{cfg.cwd}/detect/{cfg.datevshot}/\""

        if os.path.exists(f"{cfg.cwd}/vdrp/shifts/gaia/{cfg.datevshot}/all.mch"):
            cmd += f" -r \"{cfg.cwd}/vdrp/shifts/gaia\""
        elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/sdss/{cfg.datevshot}/all.mch"):
            cmd += f" -r \"{cfg.cwd}/vdrp/shifts/sdss\""
        elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/panstarrs/{cfg.datevshot}/all.mch"):
            cmd += f" -r \"{cfg.cwd}/vdrp/shifts/panstarrs\""
        else:
            print(f"[{cfg.datevshot}] Fatal: cannot find *.dithall file for this shot")
            return -1

        system_command(cfg, cmd)

        # #for convenience with legacy functions, duplicate the Shot group as a Survey group
        # shot_h5 = tables.open_file(f"{cfg.datevshot}.h5","a")
        # shot_h5.root.Shot._f_copy(newparent="root", newname="Survey", recursive=True, createparents=True)
        # shot_h5.close()


        #for convenience, softlink the shot h5 in parent directory
        #here
        shot_link_dir = os.path.join(cfg.cwd_orig, f"shots")
        if not os.path.exists(shot_link_dir):
            lock = FileLock(os.path.join(cfg.cwd_orig,Lock_mutex_fn)) #need to specify the top directory lock file
            with lock:
                Path(shot_link_dir).mkdir(parents=True, exist_ok=True)

        #assume this is the only active task working on this shot
        shot_link = os.path.join(shot_link_dir,f"{cfg.datevshot}.h5")
        shot_src = f"{cfg.cwd}/{cfg.datevshot}.h5"
        if os.path.exists(shot_link):
            system_command(cfg, f"unlink {shot_link}")

        system_command(cfg, f"ln -s {shot_src} {shot_link}")
        print(f"[{cfg.datevshot}] Created link to shot h5 file: {shot_src} -> {shot_link}")

        print(f"[{cfg.datevshot}] Done: {cfg.cwd}/{cfg.datevshot}.h5")

    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could not build hdf5 shot file.", traceback.format_exc())

    return rc


def count_amps(cfg):
    """

    :param cfg:
    :return:
    """

    bad_amps_list = []
    amps_list = []
    try:
        if os.path.exists(os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")):
            h5 = tables.open_file(os.path.join(cfg.cwd,f"{cfg.datevshot}.h5"), mode='r')
            amps_list = list(h5.root.AmpStats.read(field="multiframe").astype(str))
            bad_amps_list = list(h5.root.AmpStats.read_where("flag==0", field="multiframe").astype(str))
            cfg.num_bad_amps = len(bad_amps_list)
            cfg.num_all_amps = len(amps_list)
            h5.close()
        else:
            print(f"[{cfg.datevshot}] Could not load amps_list. shot h5 file not found.")
    except:
        print(f"[{cfg.datevshot}] Could not load amps_list",traceback.format_exc())

    return amps_list,bad_amps_list

def amp_stats(cfg,shot_h5_fqfn=None,update=True):
    """

    this needs the shot h5 file to have been completed to work
    :param cfg:
    :param shot_h5_fqfn
    :param update: if True (default) write out amp stats table. If false, recompute for internal use, but do not
                   overwrite or update existing tables
    :return:
    """

    if not update:
        print(f"[{cfg.datevshot}] Recomputing Amp Stats for internal use only. Will not update recorded tables.")

    if cfg.amp_stats_problem < 0:
        cfg.amp_stats_problem = 0

    if cfg.mf_file_status is None or len(cfg.mf_file_status.keys()) == 0:
        get_multifits_file_status(cfg)

    try:
        rc = 0

        if shot_h5_fqfn is None:
            shot_h5_fqfn = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")

        print(f"[{cfg.datevshot}] Computing amp statistics from: {shot_h5_fqfn} ... ")
        shot_dict = AmpStats.make_stats_for_shot(fqfn=shot_h5_fqfn,save=True,preload=False)

        if shot_dict is not None:
            t = AmpStats.stats_shot_dict_to_table(shot_dict)

            #t = t[t['n_lo'] >= 0] #use n_lo column to select ... the -1 values are where this failed
                                  # (e.g. usually for dithers that don't exist)

            # instead, might want to keep the n_lo == -1 ... (if they match to an expected dither)
            #      these would indicate BAD AMPS?

            sel_t = t['n_lo'] >= 0
            stat_exps = np.unique(t['expnum'][sel_t])
            #now reset and keep all with a valid exp
            sel_t = np.array([e in stat_exps for e in t['expnum']])
            t = t[sel_t]

            num_amp_stats = len(t)
            sel_failed_amp_stats = t['n_lo'] < 0
            num_failed_amp_stats = np.count_nonzero(sel_failed_amp_stats)
            unexplained_failed_amp_stats = num_failed_amp_stats
            print(f"[{cfg.datevshot}] Amp status computed for {num_amp_stats} amps, {num_failed_amp_stats} failed.")

            #are they failed amp stats okay? that is, is the multifits actually missing?
            if num_failed_amp_stats > 0:
                if cfg.mf_file_status is not None and len(cfg.mf_file_status.keys()) > 0:
                    #multiframe (s20) expnum(int64)
                    for row in t[sel_failed_amp_stats]:
                        try:
                            q_slot = row['multiframe'][10:13]
                            q_amp = row['multiframe'][18:20]
                            q_exp = row['expnum']
                            if cfg.mf_file_status[q_slot][q_exp][q_amp] != -1:
                                #the file DOES exist, so this should have generated an amp status
                                print(f"[{cfg.datevshot}] WARNING! No Amp Status, but file does exist for "
                                      f"{row['multiframe']} ({row['expnum']})")
                            else:
                                unexplained_failed_amp_stats -= 1 #makes sense, the multifits is missing
                        except:
                            print(f"[{cfg.datevshot}] WARNING! Failed mf status check for {row['multiframe']} ({row['expnum']})")
                else:
                    print(f"[{cfg.datevshot}] WARNING! Cannot confirm missing multifits for failed amp status")


            # expectation is, generally, 78 IFUs x 4 amps x # exposures = 936 for full array, 3 exposures
            #                                                            312 for full array 1 exposure
            # noting: like 095 RU might be absent, so 311 might be reasonable
            #          or earlier data might not have a full array


            if num_failed_amp_stats > 2 * len(stat_exps): #allow up to 2 (* number of exposures)
                # may be a problem
                print(f"[{cfg.datevshot}] WARNING! {num_failed_amp_stats} failures to build amp stats."
                      f" This may indicate a problem with the shot.")
                cfg.amp_stats_problem = 1

            if unexplained_failed_amp_stats != 0:
                # may be a problem
                print(f"[{cfg.datevshot}] WARNING! {unexplained_failed_amp_stats} unexplained failures to build amp stats."
                      f" This may indicate a problem with the shot.")
                cfg.amp_stats_problem = 1


            #???how much of stats_qc needs to be re-done since it is based on 3-dithers and some joint statistics???
            #several of the checks are looking for extreme variation over the dithers, which can't be done with just one dither
            try:
                t = AmpStats.stats_qc(t, extend=True,total_exp_time=cfg.total_exp_time,num_exposures=cfg.numexp)
            except: #some older versions don't have "num_exposures"
                t = AmpStats.stats_qc(t, extend=True, total_exp_time=cfg.total_exp_time)

            print(f"[{cfg.datevshot}] {np.count_nonzero(t['flag']!=1)} amps marked explicitly 'bad'")

            #for THIS (single shot reduction) we will ALSO assume the amp to be bad if n_lo < 0
            sel_bad_n_lo = t['n_lo']<0
            if np.count_nonzero(sel_bad_n_lo) > 0:
                t['flag'][sel_bad_n_lo] = 0 #0 is bad
                print(f"[{cfg.datevshot}] {np.count_nonzero(sel_bad_n_lo)} additional amps marked bad for failure"
                      f" to compute stats.")


            if update:
                t.write(f"{cfg.datevshot}_ampstats.fits",format="fits",overwrite=True)
                t.write(f"{cfg.datevshot}_ampstats.tab", format="ascii",overwrite=True)

                #always creat the bad amps file, even if none trigger
                with open(f"{cfg.datevshot}_badamps.txt","w") as f:
                    for row in t[t['flag']!=1]: #1 == Good, 0 = Bad
                        f.write(f"{row['multiframe']} exp{str(row['expnum']).zfill(2)} \n")

            cfg.num_bad_amps = np.count_nonzero(t['flag']!=1)
            cfg.num_all_amps = len(t)

            #assuming these are post-HETDEX, go ahead and put this in the shot.h5 file
            #NOTICE: this is not done for the original HETDEX shots
            #needs the actual h5 file

            #NOTICE: the "flag" key DOES NOT EXIST here ... it is added to table t above, but not to the shot_dict
            if update:
                print(f"[{cfg.datevshot}] Adding AmpStats table to  {shot_h5_fqfn} ...")
                h5 = tables.open_file(shot_h5_fqfn,mode="a")

                try:
                    #need to update the shot_dict with the results from stats_qc
                    AmpStats.stats_update_shot(h5,shot_dict=None,shot_dict_tab=t)
                except:
                    print(f"Unable to update amp stats. Exception:",traceback.format_exc())

                h5.close()

            del t

        else:
            print(f"[{cfg.datevshot}] FAIL. Could not compute amp stats.")
            cfg.amp_stats_problem = 2
            rc = -1


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could not produce amp statistics.", traceback.format_exc())
        cfg.amp_stats_problem = 3

    return rc

def add_fiber_index(cfg,shot_h5_fqfn=None):
    """

    :param cfg:
    :return:
    """

    try:
        rc = 0
        if shot_h5_fqfn is None:
            shot_h5_fqfn = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")

        #this might be an update of masking to existing shot, .... if the FiberIndex already exists, we
        #can skip this part

        h5 = tables.open_file(shot_h5_fqfn,mode='r')
        if h5.__contains__("/FiberIndex"):
            print(f"[{cfg.datevshot}] FiberIndex already exists in {shot_h5_fqfn} ... ")
            h5.close()
            return 0
        h5.close() #must close the read handle, since the create_fiber_index_hdf5 will open as append

        print(f"[{cfg.datevshot}] Adding fiber index  to: {shot_h5_fqfn} ... ")

        cmd = f"python3 {hetdex_api_path}/h5tools/create_fiber_index_hdf5.py"
        cmd += f" --shot_h5  {shot_h5_fqfn}"

        system_command(cfg, cmd)

        #test:
        try:
            h5 = tables.open_file(shot_h5_fqfn,mode="r")
            if h5.__contains__("/FiberIndex"):
                #success
                print(f"[{cfg.datevshot}] Pass. Successfully added FiberIndex.")
            else:
                print(f"[{cfg.datevshot}] FAIL! Failed to add FiberIndex.")
                rc = -1
            h5.close()

        except:
            rc = -1
            print(f"[{cfg.datevshot}] Could add fiber index to shot h5.", traceback.format_exc())
    except:
        rc = -1
        print(f"[{cfg.datevshot}] Could add fiber index to shot h5.", traceback.format_exc())

    return rc

def add_fiber_mask(cfg,shot_h5_fqfn=None):
    """

    :param cfg:
    :return:
    """

    try:
        rc = 0
        if shot_h5_fqfn is None:
            shot_h5_fqfn = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")

        print(f"[{cfg.datevshot}] Adding fiber masking (CalfibDQ) to: {shot_h5_fqfn} ... ")

        cmd = f"python3 {hetdex_api_path}/h5tools/create_fiber_mask_hdf5.py"
        cmd += f" --shot_h5  {shot_h5_fqfn}"

        system_command(cfg, cmd)

        #test:
        try:
            h5 = tables.open_file(shot_h5_fqfn,mode="r")
            if h5.__contains__("/CalfibDQ"):
                #success
                print(f"[{cfg.datevshot}] Pass. Successfully added CalfibDQ.")
            else:
                print(f"[{cfg.datevshot}] FAIL! Failed to add CalfibDQ.")
                rc = -1
            h5.close()

        except:
            rc = -1
            print(f"[{cfg.datevshot}] Could add fiber masking to shot h5.", traceback.format_exc())
    except:
        rc = -1
        print(f"[{cfg.datevshot}] Could add fiber masking to shot h5.", traceback.format_exc())

    return rc


def make_avg_mf_row(cfg):
    """

    make an "Average" row for an amp multiframe

    take median down each column for each ifu+amp+exp
    then take the mean of those medians

    result is a 1032 array of floats as an "average" row

    :param cfg:
    :return:
    """

    img = "clean_image"
    row_list = []
    avg_row = None

    use_shot_h5 = False
    shot_h5_fqfn = os.path.join(cfg.cwd, f"{cfg.datevshot}.h5")

    #for debug
    if not os.path.exists(shot_h5_fqfn):
        shot_h5_fqfn = os.path.join(cfg.cwd, f"sci{cfg.datevshot}/{cfg.datevshot}.h5")

    if os.path.exists(shot_h5_fqfn):
        try:
            h5 = tables.open_file(shot_h5_fqfn)
            data_array = h5.root.Data.Images.read(field=img)
            # there can be duplicates
            data_mf = h5.root.Data.Images.read(field="multiframe")
            data_exp = h5.root.Data.Images.read(field="expnum")

            data_mfe = [m.decode() + "_" + str(e) for m, e in zip(data_mf, data_exp)]

            umf, imf, cmf = np.unique(data_mfe, return_index=True,return_counts=True)
            data_array = data_array[imf]

            for data in data_array:
                row_list.append(np.nanmedian(data, axis=0))  # median "down" each column to make an average row

            h5.close()
            use_shot_h5 = True
        except:
            use_shot_h5 = False
            try:
                h5.close()
            except:
                pass


    if not use_shot_h5:

        # get all the reduction paths
        date = cfg.datevshot[0:8]
        virus_shot = "virus0000" + cfg.datevshot[-3:]
        if red1path is not None:
            datadir = os.path.join(red1path, f"{date}/virus/{virus_shot}")  # /{exp}/virus")
        else:
            datadir = os.path.join(cfg.cwd_orig, f"reductions/{date}/virus/{virus_shot}")  # /{exp}/virus")

        exps = sorted(glob.glob(f"{datadir}/exp*"))
        expdirs = [d + "/virus" for d in exps]

        mfs_in_exp = []
        for expdir in expdirs:
            mfs_in_exp += list(sorted(glob.glob(f"{expdir}/multi_*.fits")))


        row_list = []
        avg_row = None
        for mf_fn in mfs_in_exp:
            try:
                with fits.open(mf_fn) as hdu:
                    data = hdu[img].data
                    row_list.append(np.nanmedian(data,axis=0)) #median "down" each column to make an average row
            except:
                pass


    #now take the mean of those rows
    if len(row_list) > 30: #should normally be 312 per exposure, but early data might have only 9 IFUS or 36 amps
        avg_row = np.nanmean(row_list,axis=0)

    #for safety:
    avg_row[avg_row == 0] = 1e-6 #small, nonzero count

    return avg_row




def shot_analyisis(cfg,ratio=False):
    """

    this needs the shot h5 file to have been completed to work
    :param cfg:
    :return:
    """

    try:
        print(f"[{cfg.datevshot}] Making basic IFU analysis images ...")

        # mf_avg_row = None
        # if residual:
        #     if cfg.mf_clean_image_avg_row is None:
        #         cfg.mf_clean_image_avg_row = make_avg_mf_row(cfg)
        #     mf_avg_row = cfg.mf_clean_image_avg_row
        # else:
        #     mf_avg_row = None

        if ratio:
            print(f"[{cfg.datevshot}] NOTICE: shot_analysis with ratio = True, not corrently supported. "
                  f"'processed' image not available in shot h5 file.")
            ratio = False

        rc = 0
        shot_h5_fqfn = os.path.join(cfg.cwd, f"{cfg.datevshot}.h5")
        h5 = tables.open_file(shot_h5_fqfn)

        #make the output location
        analysis_dir = os.path.join(cfg.cwd,"analysis/")
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(analysis_dir)

        ######################################################
        # matched is already present in /matched_pngs directory
        # but can make a link
        ######################################################
        print(f"[{cfg.datevshot}] Linking match_png")
        system_command(cfg, "ln -s ../match_pngs/match*.png .")


        ######################################################
        # coadds ... just put at the top of /analysis
        ######################################################
        print(f"[{cfg.datevshot}] Making coadd image(s): ... ")

        try:
            plt.close('all')
            plt.figure(figsize=(8, 8))
            plt.title(f"{cfg.datevshot} x1")
            plt.imshow(h5.root.Astrometry.CoaddImages.png_exp01.read())  # ,origin="upper")
            plt.tight_layout()
            plt.savefig(f"coadd_{cfg.datevshot}_exp01.png", dpi=96)
        except:
            pass

        try:
            plt.close('all')
            plt.figure(figsize=(8, 8))
            plt.title(f"{cfg.datevshot} x2")
            plt.imshow(h5.root.Astrometry.CoaddImages.png_exp02.read())  # ,origin="upper")
            plt.tight_layout()
            plt.savefig(f"coadd_{cfg.datevshot}_exp02.png", dpi=96)
        except:
            pass

        try:
            plt.close('all')
            plt.figure(figsize=(8, 8))
            plt.title(f"{cfg.datevshot} x3")
            plt.imshow(h5.root.Astrometry.CoaddImages.png_exp03.read())  # ,origin="upper")
            plt.tight_layout()
            plt.savefig(f"coadd_{cfg.datevshot}_exp03.png", dpi=96)
        except:
            pass

        ##################################################
        # 4amp x 3dither (normally) pngs for each IFU
        # there are many, so make a subdir
        ##################################################
        if cfg.made_amp_images:
            print(f"[{cfg.datevshot}] Diagnostic IFU+amp images already created. Will not re-run here.")
        else:
            analysis_dir = os.path.join(cfg.cwd,"analysis/ifus/")
            Path(analysis_dir).mkdir(parents=True, exist_ok=True)
            os.chdir(analysis_dir)
            mfs_in_shot = np.unique(h5.root.Data.FiberIndex.read(field="multiframe"))
            ifus_in_shot = np.unique([mfs[0:-3] for mfs in mfs_in_shot])
            slotids = [ifu[10:13] for ifu in ifus_in_shot]
            #put these in order of slotid
            slotids, ifus_in_shot = zip(*sorted(zip(slotids, ifus_in_shot)))

            img = "clean_image"
            #cmap = "gray"
            cmap = "Spectral" #try a diverging colormap to help extremes stand out
            #use this to force the colormap to be centered at a value of 0

            # !!! NOTICE ... see also make_amp_images() for the almost identical code which can be called instead !!!!!

            #todo: this should be scaled based on the avg sky? and or time? ... do we know that yet?
            # typical 360 sec avg sky is in the 100-ish range (say 100-300 would not be uncommon)
            #  depends on moon, etc
            # HOWEVER, while the average seems to increase with sky and the stddev broadens maybe with sqrt(time)
            #   we still get negative values, so, maybe we do nothing??
            # NOTICE: the negative end does not really move, only the positive end (vmax)

            applied_vmin = int(DIAG_AMP_IMG_VMIN_VMAX[0])
            applied_vmax = int(DIAG_AMP_IMG_VMIN_VMAX[1])
            if ratio is False:
                try:
                    #vmax_scale = min(1.0,np.sqrt(cfg.total_exp_time / (360. * cfg.numexp)))
                    #no ... more linear with time, rather than sqrt, but the avg sky seems to plays a part too (stretching
                    # out more than just the time) ... maybe linear time stretch + avg_sky stretch??
                    # the "zero" (average peak) seems consistently just about 0 cts, regardless, with a positive skew
                    #  as you would expect
                    #sky_stretch = cfg.avg_sky / 300.0 #where 300 is a typical-ish HETDEX sky
                    vmax_scale = min(1.0, cfg.total_exp_time / 360.0 /cfg.numexp)
                    if vmax_scale == 0:
                        vmax_scale = 1.0
                except:
                    vmax_scale = 1.0

                #vmin_vmax_shift = 0.0 #we do not want to shift ... only scale (stretch) the vmax (positive) side
                #not sure we want to shift ... maybe just scale with time


                #print(f"[{cfg.datevshot}] Shifting IFU diagnostic plot scaling by +{vmax_scale:0.1f}")
                applied_vmin = int(DIAG_AMP_IMG_VMIN_VMAX[0])
                applied_vmax = int(DIAG_AMP_IMG_VMIN_VMAX[1] * vmax_scale)
                cmap_norm = TwoSlopeNorm(vmin=applied_vmin, vcenter=0, vmax=applied_vmax)
            else:
                vmax_scale = 1.0
                applied_vmin = -1
                applied_vmax = 1
                cmap_norm = TwoSlopeNorm(vmin=applied_vmin, vcenter=0, vmax=applied_vmax)


            #just assume 3 dithers ... if they do not exist, they will be blank
            #there are a few that have 4 or more exposures and we will just ignore that
            for ii, mf_base in enumerate(ifus_in_shot):
                if ratio:
                    print(f"[{cfg.datevshot}] Making ratio IFU analysis images: {mf_base.decode()}")
                else:
                    print(f"[{cfg.datevshot}] Making basic IFU analysis images: {mf_base.decode()}")

                plt.close('all')
                fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(9, 12))
                #plot_config = list(np.arange(431, 443, 1))
                if ratio is False:
                    fig.suptitle(f"{cfg.datevshot} {mf_base} {img} counts, cmap scale: "
                                 f"({applied_vmin},"
                                 f" {applied_vmax})")
                else:
                    fig.suptitle(f"{cfg.datevshot} {mf_base} {img}/processed, cmap scale: "
                                 f"({applied_vmin},"
                                 f" {applied_vmax})")

                for ai, amp in enumerate([b'_RU', b'_RL', b'_LL', b'_LU']):
                    ei = 0  # exposure index
                    # print(ai,ei,amp)
                    mf = mf_base + amp
                    data = h5.root.Data.Images.read_where("multiframe==mf", field=img)
                    if ratio:
                        proc = h5.root.Data.Images.read_where("multiframe==mf", field="processed")
                        data = data / (proc-data)
                    exp = h5.root.Data.Images.read_where("multiframe==mf", field='expnum')

                    # ax = plt.subplot(plot_config[ci])
                    ax = axes[ai, ei]
                    ei += 1
                    sel = exp == 1
                    try:
                        if np.count_nonzero(sel) == 1:
                            ax.set_title(f"{amp.decode()} x1")
                            vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                            #cmap_norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                            if amp != b'_LU':
                                ax.set_xticks([])
                            else:
                                ax.set_xticks([0,200,400,600,800,1000])
                            #always set yticks on 1st column
                            ax.set_yticks([0, 200, 400, 600, 800, 1000])
                            # !!! Must use EITHER norm or vmin, vmax ... cannot do both
                            #ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax,origin="lower")
                            ax.imshow(data[sel][0], cmap=cmap, norm=cmap_norm, origin="lower")
                        else:
                            ax.set_xticks([])
                            ax.set_yticks([])
                            #not uncommon ... could be just a single exposure by user selection, and not exp 1
                            #print(f"[{cfg.datevshot}] No data found for {mf_base.decode()} exp 1")

                    except:
                        ax.set_xticks([])
                        ax.set_yticks([])
                        print(f"[{cfg.datevshot}] Exception in shot_analysis on {mf_base.decode()}", traceback.format_exc())

                    sel = exp == 2
                    ax = axes[ai, ei]
                    ei += 1
                    try:
                        if np.count_nonzero(sel) == 1:

                            ax.set_title(f"{amp.decode()} x2")
                            vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                            #cmap_norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                            # ax.yaxis.label.set_visible(False)

                            if amp != b'_LU':
                                ax.set_xticks([])
                            else:
                                ax.set_xticks([0,200,400,600,800,1000])

                            # always unset yticks on not 1st column
                            ax.set_yticks([])
                            # !!! Must use EITHER norm or vmin, vmax ... cannot do both
                            #ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax,origin="lower")
                            ax.imshow(data[sel][0], cmap=cmap, norm=cmap_norm, origin="lower")
                        else:
                            ax.set_xticks([])
                            ax.set_yticks([])
                            #not uncommon ... could be just a single exposure for this observation or by user selection
                            #print(f"[{cfg.datevshot}] No data found for {mf_base.decode()} exp 2")

                    except:
                        ax.set_xticks([])
                        ax.set_yticks([])
                        print(f"[{cfg.datevshot}] Exception in shot_analysis on {mf_base.decode()}",
                              traceback.format_exc())

                    sel = exp == 3
                    ax = axes[ai, ei]
                    ei += 1
                    try:
                        if np.count_nonzero(sel) == 1:
                            ax.set_title(f"{amp.decode()} x3")
                            vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                            #cmap_norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                            if amp != b'_LU':
                                ax.set_xticks([])
                            else:
                                ax.set_xticks([0,200,400,600,800,1000])

                            # always unset yticks on not 1st column
                            ax.set_yticks([])
                            #!!! Must use EITHER norm or vmin, vmax ... cannot do both
                            #ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax,origin="lower")
                            ax.imshow(data[sel][0], cmap=cmap, norm=cmap_norm, origin="lower")
                        else:
                            ax.set_xticks([])
                            ax.set_yticks([])
                            #not uncommon ... could be just a single exposure for this observation or by user selection
                            #print(f"[{cfg.datevshot}] No data found for {mf_base.decode()} exp 3")

                    except:
                        ax.set_xticks([])
                        ax.set_yticks([])
                        print(f"[{cfg.datevshot}] Exception in shot_analysis on {mf_base.decode()}",
                              traceback.format_exc())

                plt.tight_layout()
                if ratio is False:
                    plt.savefig(f"i{mf_base.decode()[10:13]}_{cfg.datevshot}_{mf_base.decode()}.png",
                                dpi=DIAG_AMP_IMG_DPI)
                else:
                    plt.savefig(f"i{mf_base.decode()[10:13]}_{cfg.datevshot}_{mf_base.decode()}_ratio.png",
                                dpi=DIAG_AMP_IMG_DPI)

                cfg.made_amp_images = True
    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could complete basic shot analysis.",traceback.format_exc())

    try:
        h5.close()
    except:
        pass

    return rc

def get_multifits_file_status(cfg):
    """
    check the reduction diretory and make a dictionary with key = multifits name
    and value = -1 (missing) 0 = okay, 1 = damaged/warning
    for the multifits file under the reductions directory

    :param cfg:
    :return:
    """
    d = {}
    try:
        print(f"[{cfg.datevshot}] Getting multi*fits file status ...")

        #slotids, show 010 to 109, noting only 78 of them actually exist
        dkeys = {
           # '010',
           # '011',
           # '012',
            '013',
            '014',
            '015',
            '016',
            #'017',
            #'018',
            #'019',

            #'020',
            '021',
            '022',
            '023',
            '024',
            '025',
            '026',
            '027',
            '028',
            #'029',

            '030',
            '031',
            '032',
            '033',
            '034',
            '035',
            '036',
            '037',
            '038',
            '039',

            '040',
            '041',
            '042',
            '043',
            '044',
            '045',
            '046',
            '047',
            '048',
            '049',

            '050',
            '051',
            '052',
            '053',
            #'054',
            #'055',
            #'056',
            '057',
            '058',
            '059',

            '060',
            '061',
            '062',
            '063',
            #'064',
            #'065',
            #'066',
            '067',
            '068',
            '069',

            '070',
            '071',
            '072',
            '073',
            '074',
            '075',
            '076',
            '077',
            '078',
            '079',

            '080',
            '081',
            '082',
            '083',
            '084',
            '085',
            '086',
            '087',
            '088',
            '089',

            #'090',
            '091',
            '092',
            '093',
            '094',
            '095',
            '096',
            '097',
            '098',
            #'099',

            # '100',
            # '101',
            # '102',
            '103',
            '104',
            '105',
            '106',
            # '107',
            # '108',
            # '109',
        }

        #get all the reduction paths
        date = cfg.datevshot[0:8]
        virus_shot = "virus0000" + cfg.datevshot[-3:]
        if red1path is not None:
            datadir = os.path.join(red1path, f"{date}/virus/{virus_shot}") #/{exp}/virus")
        else:
            datadir = os.path.join(cfg.cwd_orig, f"reductions/{date}/virus/{virus_shot}") #/{exp}/virus")

        exps = sorted(glob.glob(f"{datadir}/exp*"))
        expdirs = [d + "/virus" for d in exps]
        exps = [int(os.path.basename(d)[-2:]) for d in exps] #now an integer 1, 2, 3 ...

        d = {key: {} for key in dkeys}
        for exp, expdir in zip(exps,expdirs):
            for k in dkeys:
                d[k][exp]={'LL':-1,'LU':-1,'RL':-1,'RU':-1} #initialize to missing (-1)

            mfs_in_exp = [os.path.basename(x) for x in glob.glob(f"{expdir}/multi_*.fits")]
            for mf in mfs_in_exp: #i.e. multi_421_069_006_LU.fits
                slot = mf[10:13]
                amp = mf[18:20]
                d[slot][exp][amp] = 0

        #any remaing -1 are missing
        cfg.mf_file_status = d
    except:
        print(f"[{cfg.datevshot}] Could not get multi*fits file status",traceback.format_exc())

    return d

def make_amp_images(cfg,ratio=False):
    """

    originally shot_analysis ... and is based on that but intended only for the
    color-coded amp images and to be run at the end of step1
    :param cfg:
    :param residual:
    :return:
    """

    try:
        print(f"[{cfg.datevshot}] Making color-coded IFU+amp images ...")

        # mf_avg_row = None
        # if residual:
        #     if cfg.mf_clean_image_avg_row is None:
        #         cfg.mf_clean_image_avg_row = make_avg_mf_row(cfg)
        #     mf_avg_row = cfg.mf_clean_image_avg_row
        # else:
        #     mf_avg_row = None

        rc = 0

        #make the output location
        analysis_dir = os.path.join(cfg.cwd,"analysis/")
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(analysis_dir)

        analysis_dir = os.path.join(cfg.cwd, "analysis/ifus/")
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(analysis_dir)


        #get all the reduction paths
        date = cfg.datevshot[0:8]
        virus_shot = "virus0000" + cfg.datevshot[-3:]
        if red1path is not None:
            datadir = os.path.join(red1path, f"{date}/virus/{virus_shot}") #/{exp}/virus")
        else:
            datadir = os.path.join(cfg.cwd_orig, f"reductions/{date}/virus/{virus_shot}") #/{exp}/virus")

        exps = sorted(glob.glob(f"{datadir}/exp*"))
        expdirs = [d + "/virus" for d in exps]
        exps = [int(os.path.basename(d)[-2:]) for d in exps] #now an integer 1, 2, 3 ...

        if len(expdirs) == 1: #this is common
            expdirs.append("")
            expdirs.append("")
        elif len(expdirs) == 2:
            expdirs.append("")


        #could assume the same mfs in each directory, though there could be some weird failure

        mfs_in_exp = []
        for expdir in expdirs:
            mfs_in_exp += list(sorted(glob.glob(f"{expdir}/multi_*.fits")))

        #strip off the _RL.fits
        ifus_in_exp = np.unique([os.path.basename(mfs)[0:-8] for mfs in mfs_in_exp])
        slotids = [ifu[10:13] for ifu in ifus_in_exp]
        #put these in order of slotid
        slotids, ifus_in_exp = zip(*sorted(zip(slotids, ifus_in_exp)))

        img = "clean_image"
        #cmap = "gray"
        cmap = "Spectral" #try a diverging colormap to help extremes stand out
        # todo: this should be scaled based on the avg sky? and or time? ... do we know that yet?
        # typical 360 sec avg sky is in the 100-ish range (say 100-300 would not be uncommon)
        #  depends on moon, etc
        # HOWEVER, while the average seems to increase with sky and the stddev broadens maybe with sqrt(time)
        #   we still get negative values, so, maybe we do nothing??
        # NOTICE: the negative end does not really move, only the positive end (vmax)

        #!!! NOTICE ... see also shot_analysis() for the almost identical code which can be called instead !!!!!
        applied_vmin = int(DIAG_AMP_IMG_VMIN_VMAX[0])
        applied_vmax = int(DIAG_AMP_IMG_VMIN_VMAX[1])
        if ratio is False:
            try:
                # vmax_scale = min(1.0,np.sqrt(cfg.total_exp_time / (360. * cfg.numexp)))
                # no ... more linear with time, rather than sqrt, but the avg sky seems to plays a part too (stretching
                # out more than just the time) ... maybe linear time stretch + avg_sky stretch??
                # the "zero" (average peak) seems consistently just about 0 cts, regardless, with a positive skew
                #  as you would expect
                # sky_stretch = cfg.avg_sky / 300.0 #where 300 is a typical-ish HETDEX sky
                vmax_scale = min(1.0, cfg.total_exp_time / 360.0 / cfg.numexp)
                if vmax_scale == 0:
                    vmax_scale = 1.0
            except:
                vmax_scale = 1.0

            # vmin_vmax_shift = 0.0 #we do not want to shift ... only scale (stretch) the vmax (positive) side
            # not sure we want to shift ... maybe just scale with time

            # print(f"[{cfg.datevshot}] Shifting IFU diagnostic plot scaling by +{vmax_scale:0.1f}")
            # print(f"*** DEBUG *** vmin = {int(DIAG_AMP_IMG_VMIN_VMAX[0])}, "
            #       f"vmax = {int(DIAG_AMP_IMG_VMIN_VMAX[1] * vmax_scale)}, vmax_scale = {vmax_scale}, "
            #       f"exptime = {cfg.total_exp_time}, expnum = {cfg.numexp}")
            applied_vmin = int(DIAG_AMP_IMG_VMIN_VMAX[0])
            applied_vmax = int(DIAG_AMP_IMG_VMIN_VMAX[1] * vmax_scale)
            cmap_norm = TwoSlopeNorm(vmin=applied_vmin, vcenter=0,vmax=applied_vmax)
        else:
            vmax_scale = 1.0
            applied_vmin = -1
            applied_vmax = 1
            cmap_norm = TwoSlopeNorm(vmin=applied_vmin, vcenter=0,vmax=applied_vmax)

        #just assume 3 dithers ... if they do not exist, they will be blank
        #there are a few that have 4 or more exposures and we will just ignore that
        for ii, mf_base in enumerate(ifus_in_exp):

            handy_im = None
            handy_ax = None

            if ratio:
                print(f"[{cfg.datevshot}] Making ratio IFU analysis images: {mf_base}")
            else:
                print(f"[{cfg.datevshot}] Making basic IFU analysis images: {mf_base}")
            plt.close('all')
            fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(9, 12))
           # plot_config = list(np.arange(431, 443, 1))
            if ratio:
                fig.suptitle(f"{cfg.datevshot} {mf_base} {img}/(processed-{img}), cmap scale: "
                             f"({applied_vmin},"
                             f" {applied_vmax})\n\n")
            else:
                fig.suptitle(f"{cfg.datevshot} {mf_base} {img} counts, cmap scale: "
                             f"({applied_vmin},"
                             f" {applied_vmax})\n\n")


            for ei, expdir in enumerate(expdirs):

                ##################################################
                # 4amp x 3dither (normally) pngs for each IFU
                # there are many, so make a subdir
                ##################################################

                if ei > 2:
                    break #we're done for the images, limited to only 3 exp max

                try:
                    exp = int(expdir.split("/exp")[1][:2])
                except:
                    exp = ei + 1


                # images should be under the reduction directory
                for ai, amp in enumerate(['_RU', '_RL', '_LL', '_LU']):
                    # print(ai,ei,amp)
                    mf_fn = mf_base + amp + ".fits"

                    try:
                        with fits.open(os.path.join(expdir,mf_fn)) as hdu:
                            data = hdu[img].data
                            if ratio:
                                #data = (data - mf_avg_row) / mf_avg_row
                                proc = hdu["processed"].data
                                #proc[proc==0] = np.inf
                                data = data / (proc-data)
                            ax = axes[ai, ei]

                            try:
                                ax.set_title(f"{amp} x{exps[ei]}")
                                #vmin, vmax = Utils.get_vrange(data, contrast=0.25)
                                #cmap_norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
                                if amp == '_LU':
                                    ax.set_xticks([0, 200, 400, 600, 800, 1000])
                                else:
                                    ax.set_xticks([])

                                if ei == 0:
                                    #always set yticks on 1st column
                                    ax.set_yticks([0, 200, 400, 600, 800, 1000])
                                else:
                                    # always unset yticks on not 1st column
                                    ax.set_yticks([])
                                # !!! Must use EITHER norm or vmin, vmax ... cannot do both
                                #ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax,origin="lower")
                                handy_im = ax.imshow(data, cmap=cmap, norm=cmap_norm, origin="lower")
                                handy_ax = ax

                            except:
                                ax.set_xticks([])
                                ax.set_yticks([])
                                print(f"[{cfg.datevshot}] Exception in make_amp_images on {mf_base}", traceback.format_exc())
                    except:
                        if os.path.exists(os.path.join(expdir, mf_fn)):
                            print(f"[{cfg.datevshot}] Exception in make_amp_images on {mf_base}",
                                  traceback.format_exc())
                        else:
                            if len(expdir) > 1:
                                #the multifits does not exist, probably a pre-defined bad amp like 095RU
                                #if expdir is an empty string, then it is just there to clean up the
                                #  image frames and was not to be printed
                                print(f"[{cfg.datevshot}] omitting {mf_base}{amp} x{exp}. Does not exist.")
                            try:
                                ax = axes[ai, ei]
                                if amp == '_LU':
                                    ax.set_xticks([0, 200, 400, 600, 800, 1000])
                                else:
                                    ax.set_xticks([])

                                if ei == 0:
                                    #always set yticks on 1st column
                                    ax.set_yticks([0, 200, 400, 600, 800, 1000])
                                else:
                                    # always unset yticks on not 1st column
                                    ax.set_yticks([])
                            except:
                                pass


            #fig.colorbar(handy_im, ax=axes, location='bottom', cmap=cmap, norm=cmap_norm, pad=-0.45)#,
            # axins = inset_axes(
            #     axes,
            #     height="5%",  # Width of the colorbar (relative to the bbox)
            #     width="100%",  # Height of the colorbar (relative to the bbox)
            #     loc="lower left",  # Alignment point inside the bbox
            #     bbox_to_anchor=(1.05, 0., 1, 1),  # (x, y, width, height) relative to transAxes
            #     bbox_transform=handy_ax.transAxes,  # Use normalized axes coordinates (0 to 1)
            #     borderpad=0  # Remove padding around the inset axes
            # )
            #
            # fig.colorbar(handy_im, cax=axins, location='bottom', anchor=(0.5,1.0),panchor=(0.5,1.0),
            #              fraction=0.1, cmap=cmap, norm=cmap_norm,pad=0.15)

            cax = fig.add_axes([0.15, 0.95, 0.7, 0.015])
            fig.colorbar(handy_im, cax=cax, orientation='horizontal',cmap=cmap, norm=cmap_norm, pad=0.0)

            #anchor=(0.0,0.0),, pad=0.15)  # ,

            plt.tight_layout()
            if ratio:
                plt.savefig(f"i{mf_base[10:13]}_{cfg.datevshot}_{mf_base}_ratio.png", dpi=DIAG_AMP_IMG_DPI)
            else:
                plt.savefig(f"i{mf_base[10:13]}_{cfg.datevshot}_{mf_base}.png", dpi=DIAG_AMP_IMG_DPI)

    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could not complete amp images.",traceback.format_exc())


    if rc >= 0:
        cfg.made_amp_images = True #this is checked later in shot_analysis to not re-do the images IF they
                                   #were done in THIS call (note: that if this is a --resume, this will
                                   #not be set and they will be re-made, and that is on purpose
                                   #as something could have changed, including the analysis code
        #print("*** DEBUG *** forcing cfg.made_amp_images to False so the shot_analysis version will run.")
        #cfg.made_amp_images = False
    return rc


def check_detid_counter(cfg,cont,nominal_max):
    """

    :param cfg:
    :param cont: if True, check the continuum h5, else check the line h5 file
    :return: 0 if okay (has not exceeded the limit) 1 if exceeded limit, -1 is error
    """

    try:
        if cont:
            which = "cont"
        else:
            which = "line"

        if os.path.exists(f"{cfg.datevshot}_{which}.h5"):
            h5 = tables.open_file(f"{cfg.datevshot}_{which}.h5")
            ids = h5.root.Detections.read(field="detectid")
            if ids is not None and len(ids) > 0:
                maxid = np.max(ids)

                if maxid > nominal_max:
                    rc = len(ids)
                else:
                    rc = 0
            else:
                print(f"[{cfg.datevshot}] Warning! No detections found in {cfg.datevshot}_{which}.h5")
                rc = -1 #it is not reasonable that there should not be any detections
                        #this is not fatal, but it will print an error
        else:
            print(f"[{cfg.datevshot}] Warning! File does not exist: {cfg.datevshot}_{which}.h5")
            return -1 #the file should exist
    except:
        rc = -1
        print(f"[{cfg.datevshot}] Exception in check_detid_counter()",traceback.format_exc())

    try:
        h5.close()
    except:
        pass

    return rc

def build_detection_hdf5(cfg):
    """

    line detects and continuum

    build a catalog astropy table (fits format) for the line and continuum detections

    :param cfg:
    :return:
    """


    try:
        rc = 0
        id_ct = -1
        id_base = 5 #5 character counter (was 1e5 now the string version with 0's id_base + 1 long,

        while id_ct != 0:
            #strip off the first 2 digits of the year and the 'v'
            #leaves last 5 characters for a counter for the detectids in the shot
            #note: 9xxxx is reserved for continuum sources
            #detectid_base = np.int64(cfg.datevshot[2:].replace('v', '', )) * np.int64(1e5) # + 1000000
            detectid_base = np.int64(cfg.datevshot[2:].replace('v', '', )) * np.int64(str(1).ljust(id_base+1,'0'))  # + 1000000
            print(f"[{cfg.datevshot}] Building detections hdf5 ... ")

            os.chdir(os.path.join(cfg.cwd))

            cmd = f"python3 {hetdex_api_path}/h5tools/create_detect_hdf5.py"
            #cmd = f"python3 {cfg.cwd}/../create_detect_hdf5.py"
            #cmd += f" --survey NA"
            cmd += f" --detectid_base {detectid_base}"
            cmd += f" --date {cfg.datevshot[0:8]}"
            cmd += f" --observation \"{cfg.datevshot[-3:]}\""
            cmd += f" -of \"{cfg.datevshot}_line.h5\""
            cmd += f" --detect_path \"{cfg.cwd}/alldet/detect_out\""
            if cfg.linedet_filter >= 1.1:
                cmd += f" --sn_min {cfg.linedet_filter:0.1f}"
            elif cfg.snr_rescale > 1.0:
                #NOTE: the HETDEX_API minimum value is configured at 4.5 NOT 4.8
                # which is the defulat for HETDEX_API_SNR_Thresh
                cmd += f" --sn_min {HETDEX_API_SNR_Thresh * cfg.snr_rescale:0.1f}"
                print(f"[{cfg.datevshot}] Rescaled minimum SNR for inclusion from "
                      f"{HETDEX_API_SNR_Thresh:0.1f} to {HETDEX_API_SNR_Thresh * cfg.snr_rescale:0.1f}.")
            #else, using HETDEX_API default if or 1

            system_command(cfg, cmd)


            #sanity check the counter ... IF we exceed 10,000 detections, the counter fails
            id_ct = check_detid_counter(cfg,cont=False,nominal_max= detectid_base + np.int64(str(8).ljust(id_base ,'9')))
            if id_ct > 0:
                #this is a problem and we need to re-run with a different base
                print(f"[{cfg.datevshot}] More than expected ({id_ct}) number of line detections. Must reset counter and rebuild.")
                id_base = len(str(id_ct)) + 1
            elif id_ct < 0:
                #this is an error
                print(f"[{cfg.datevshot}] Error ({id_ct}) checking number of line detections. ")
                id_ct = 0


        id_ct = -1
        id_base = 5 #5 character counter (was 1e5 now the string version with 0's id_base + 1 long,

        while id_ct != 0:

            #detectid_base = np.int64(cfg.datevshot.replace('v', '', )) * np.int64(1e7) + 9000000
            #detectid_base = np.int64(cfg.datevshot[2:].replace('v', '', )) * np.int64(1e5) + 90000
            detectid_base = np.int64(cfg.datevshot[2:].replace('v', '', )) * np.int64(str(1).ljust(id_base + 1, '0')) + \
                            np.int64(str(9).ljust(id_base , '0'))
            cmd = f"python3 {hetdex_api_path}/h5tools/create_cont_hdf5.py"
            #cmd = f"python3 {cfg.cwd}/../create_cont_hdf5.py"
            cmd += f" --detectid_base {detectid_base}"
            cmd += f" --date {cfg.datevshot[0:8]}"
            cmd += f" --observation \"{cfg.datevshot[-3:]}\""
            cmd += f" -of \"{cfg.datevshot}_cont.h5\""
            cmd += f" --detect_path \"{cfg.cwd}/cs/spec\""

            system_command(cfg, cmd)

            #sanity check the counter ... IF we exceed 10,000 detections, the counter fails
            id_ct = check_detid_counter(cfg,cont=True,nominal_max=detectid_base + np.int64(str(9).ljust(id_base ,'9')))
            if id_ct > 0:
                #this is a problem and we need to re-run with a different base
                print(f"[{cfg.datevshot}] More than expected ({id_ct}) number of continuum detections. Must reset counter and rebuild.")
                id_base = len(str(id_ct)) + 1
            elif id_ct < 0:
                #this is an error
                print(f"[{cfg.datevshot}] Error ({id_ct}) checking number of continuum detections. ")
                id_ct = 0


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Exception building detections hdf5 files",traceback.format_exc())

    return rc


def diagnose(cfg):
    """

    prep and run Diagnose

    all continuum sources
    all line sources that are brighter than gmag 23 ?

    :param cfg:
    :return:
    """

    try:

        os.chdir(os.path.join(cfg.cwd))
        rc = 0
        num_spectra = 0

        print(f"[{cfg.datevshot}] Building Diagnose input (emission lines)  ...")
        name = f"{cfg.datevshot}_line_sourcecat.tab"
        line_tab = Table.read(name,format="ascii")

        #looks like MW dust correction is needed here !!!
        #it may be sufficient just use the corretion for the shot? or should we use each one individually?
        ra = None
        dec = None
        try:
            if cfg.shot_ra is not None and cfg.shot_ra > -999:
                ra = cfg.shot_ra
                dec = cfg.shot_dec
            else:
                h5 = None
                try:
                    h5 = tables.open_file(f"{cfg.datevshot}.h5", mode="r")
                    ra = h5.root.Shot.read(field='ra')[0]
                    dec = h5.root.Shot.read(field='dec')[0]
                except:
                    pass

                if h5:
                    h5.close()

            if ra is not None and dec is not None:
                dust_corr = deredden_spectra(G.CALFIB_WAVEGRID,
                                         SkyCoord(ra=ra * u.deg, dec=dec * u.deg))
                print(f"[{cfg.datevshot}] MW dust correction for Diagnose:  {dust_corr[0]:0.2f} to {dust_corr[-1]:0.2f}")
            else:
                print(f"[{cfg.datevshot}] Could not buil MW dust correction for Diagnose."
                      f". RA, Dec = ({ra},{dec})")
        except:
            print(f"[{cfg.datevshot}] Exception building MW dust correction for Diagnose."
                  f". RA, Dec = ({ra},{dec}){traceback.format_exc()}")
            dust_corr = np.ones(len(G.CALFIB_WAVEGRID))



        if len(line_tab) == 0: #there are no entries
            del line_tab
            print(f"[{cfg.datevshot}] WARNING! No line sources recored. Moving on to continuum.")
            rc = 1 #not fatal, but not a success
        else:
            line_h5 = tables.open_file(os.path.join(cfg.cwd, f"{cfg.datevshot}_line.h5"))

            #reminder: if fof clustering re-runs (the step before this one), the gmag will be wiped out
            if cfg.resume or 'gmag' not in line_tab.columns:
                print(f"[{cfg.datevshot}] Computing {len(line_tab)} line det gmags for Diagnose sub-selection ...")
                # subselect ... these do not currently have a gmag, so need to make one (since ELiXer has not run yet)
                line_tab['gmag'] = 99.0 #not-computed or failed compute value


                # pre-read all columns
                # will match to detectid
                # !!! Remember. line_tab (from the line_sourcat.tab) is downselected from line.h5
                #   ... line.h5 will have equal or MORE detections than are in line_tab
                print(f"[{cfg.datevshot}] Preloading emission line spectra ...")
                all_detid = list(line_h5.root.Spectra.read(field="detectid"))
                all_spec1d = line_h5.root.Spectra.read(field="spec1d")
                all_spec1d_err = line_h5.root.Spectra.read(field="spec1d_err")
                all_apcor = line_h5.root.Spectra.read(field="apcor")
                all_approx_gmag = np.nansum(all_spec1d,axis=1)

                #print(f"[{cfg.datevshot}] DEBUG: all_approx_gmag shape = {np.shape(all_approx_gmag)}")

                print(f"[{cfg.datevshot}] Evaluating pseudo gmags ... ")
                #for i in tqdm(range(len(line_tab))):
                for i in range(len(line_tab)):
                    d = line_tab[i]['detectid']
                    try:
                        j = all_detid.index(d) #match line_tab detectid to line.h5 detectid
                    except:
                        print(f"[{cfg.datevshot}] Warning! detectid: {d} not found in *_line.h5 file?")
                        continue

                    if all_approx_gmag[j] < 640.0: #too faint (fainter than about 23 in g)
                        #640.0 ~ sum along the flux: flux for g = 23 @ 4640AA is about 0.32e-17 erg/s/cm2/AA
                        # for 0.32 * 2AA * 1000 bins = 640.0
                        continue


                    #older method, just reading as needed ... (the newer way (above) needs more memory, but is much
                    #   faster, pre-loading all in advance)


                    # rows = line_h5.root.Spectra.read_where("detectid==d")
                    # if len(rows) != 1:
                    #     print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                    #     continue
                    # row = rows[0]
                    # gmag, gmag_unc, *_ = SU.get_best_gmag(row['spec1d'] * 1e-17, row['spec1d_err'] * 1e-17,
                    #                                      G.CALFIB_WAVEGRID)

                    #reading less but using two reads: overall time cost is about the same
                    # spec1d = line_h5.root.Spectra.read_where("detectid==d",field="spec1d")
                    # if len(spec1d) != 1:
                    #     print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                    #     continue
                    # spec1d_err = line_h5.root.Spectra.read_where("detectid==d",field="spec1d_err")
                    # gmag, gmag_unc, *_ = SU.get_best_gmag(spec1d[0]*1e-17,spec1d_err[0]*1e-17,G.CALFIB_WAVEGRID)


                    gmag, gmag_unc, *_ = SU.get_best_gmag(all_spec1d[j]*1e-17,all_spec1d_err[j]*1e-17,G.CALFIB_WAVEGRID)


                    line_tab['gmag'][i] = gmag

                #update
                line_tab.write(name, format="ascii", overwrite=True)
                print(f"[{cfg.datevshot}] Updated lines table with gmag: {os.getcwd()}/{name}")

            # now select on gmag < 23
            sel = np.array(line_tab['gmag'] > 0.0) & np.array(line_tab['gmag'] <= 23.0)
            line_tab = line_tab[sel]
            #redundant comment to one just below this section
            #print(f"[{cfg.datevshot}] Reduced Diagnose examination to {len(line_tab)} emission line detections.")


            totN = np.sum(sel)
            DETECTID = np.array(line_tab['detectid'])
            RA = np.array(line_tab['ra'])
            DEC = np.array(line_tab['dec'])
            GMAG = np.array(line_tab['gmag'])
            #WEIGHT = np.array(line_tab['apcor'])

            RMAG = np.zeros((totN,))
            IMAG = np.zeros((totN,))
            ZMAG = np.zeros((totN,))
            YMAG = np.zeros((totN,))
            SN = np.zeros((totN,))
            NAME = np.zeros((totN,))  # np.array(source_spectra['shotid'][sel])

            #re-collect the spec for those selected
            #line_tab['spec'] = np.zeros((len(line_tab),1036))
            #line_tab['spec_err'] = np.zeros((len(line_tab),1036))
            spec_2D = []
            error_2D = []
            apcor = []

            # for i in range(len(line_tab)):
            #     d = line_tab[i]['detectid']
            #
            #     rows = line_h5.root.Spectra.read_where("detectid==d")
            #     if len(rows) != 1:
            #         print(f"Diagnose preselection spectra failure for {d}")
            #         spec_2D.append(np.zeros(1036))
            #         error_2D.append(np.zeros(1036))
            #         apcor.append(np.zeros(1036))
            #         continue
            #
            #     row = rows[0]
            #     #line_tab['spec'][i] = rows['spec1d']
            #     #line_tab['spec_err'][i] = rows['spec1d_err']
            #     spec_2D.append(row['spec1d'] * dust_corr)
            #     error_2D.append(row['spec1d_err'])
            #     apcor.append(row['apcor'])

            for i in range(len(line_tab)):
                d = line_tab[i]['detectid']
                try:
                    j = all_detid.index(d)  # match line_tab detectid to line.h5 detectid
                except:
                    print(f"[{cfg.datevshot}] Warning! detectid: {d} not found in *_line.h5 file?")
                    continue

                spec_2D.append(all_spec1d[j] * dust_corr)
                error_2D.append(all_spec1d_err[j]* dust_corr)
                apcor.append(all_apcor[j])

            #SPEC = np.array(line_tab['spec'])
            #ERROR = np.array(line_tab['spec_err'])
            SPEC = np.array(spec_2D)
            ERROR = np.array(error_2D)
            WEIGHT = np.array(apcor)

            line_h5.close()

            num_spectra = len(SPEC)

            #now write out the fits file needed for Diagnose
            outname = 'diagnose_spectra_line.fits'
            T = Table([DETECTID, RA, DEC, NAME, GMAG, RMAG, IMAG, ZMAG, YMAG, SN],
                      names=['detectid', 'RA', 'Dec', 'shotid', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 'sn'])
            fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(T), fits.ImageHDU(SPEC),
                          fits.ImageHDU(ERROR), fits.ImageHDU(WEIGHT)]).writeto(outname, overwrite=True)

            #clean up
            del T
            del DETECTID
            del RA
            del DEC
            del GMAG
            del WEIGHT
            del RMAG
            del IMAG
            del ZMAG
            del YMAG
            del SN
            del NAME
            del SPEC
            del ERROR

            #call Diagnose' needs sklearn
            # looks like these count from 0, so -li is inclusive and -hi is not
            # might have changed to --low_index and --high_index  or --index_filename
            # --catalog diagnose_spectra.fits
            # --normalize     this is a switch only, multiplies the flux by some normalization ????#
            # --quick   #a switch  ... uses saved models?
            # --suffix    # probably was -s  ... writes out to "classification_XXX.fits" where XXX is the integer --suffix
            # python3 /work/05350/ecooper/stampede2/Diagnose/diagnose.py diagnose_spectra.fits -li 0 -hi 5 -s 1 -q

            #Diagnose path is: {cfg.scriptdir}/Diagnose/diagnose.py
            print(f"[{cfg.datevshot}] Running Diagnose on {num_spectra} gmag bright emission line sources ...")
            cmd = f"python3 {cfg.scriptdir}/Diagnose/diagnose.py"
            cmd += f" {outname} -li 0 -hi {num_spectra} -s 0 -q"

            system_command(cfg, cmd)
        #end bright emission line Diagnose

        print(f"[{cfg.datevshot}] Building Diagnose input (continuum)  ...")
        name = f"{cfg.datevshot}_cont_sourcecat.tab"
        cont_tab = Table.read(name, format="ascii")

        if len(cont_tab) == 0:
            del cont_tab
            print(f"[{cfg.datevshot}] WARNING! No continuum sources recored. ")
            if rc == 0: #there was not a previous problem
                rc = 1 # not fatal, but not a success
            else:
                rc = -1  #already a problem, so Diagnose fails
                print(f"[{cfg.datevshot}] WARNING! Cannot run Diagnose.")
        else:
            cont_h5 = tables.open_file(os.path.join(cfg.cwd, f"{cfg.datevshot}_cont.h5"))

            if cfg.resume or 'gmag' not in cont_tab.columns:
                #print(f"[{cfg.datevshot}] Computing gmags for Diagnose ...")
                print(f"[{cfg.datevshot}] Computing {len(cont_tab)} continuum det gmags for Diagnose sub-selection ...")
                # subselect ... these do not currently have a gmag, so need to make one (since ELiXer has not run yet)
                cont_tab['gmag'] = 99.0


                # pre-read all columns
                # will match to detectid
                # !!! Remember. cont_tab (from the cont_sourcat.tab) is downselected from cont.h5
                #   ... cont.h5 will have equal or MORE detections than are in cont_tab
                print(f"[{cfg.datevshot}] Preloading continuum spectra ...")
                all_detid = list(cont_h5.root.Spectra.read(field="detectid"))
                all_spec1d = cont_h5.root.Spectra.read(field="spec1d")
                all_spec1d_err = cont_h5.root.Spectra.read(field="spec1d_err")
                # No. If already marked as a continuum source, use Diagnose and need the real gmag
                #all_approx_gmag = np.nansum(all_spec1d,axis=1)

                print(f"[{cfg.datevshot}] Evaluating pseudo gmags ... ")
                for i in range(len(cont_tab)):
                    d = cont_tab[i]['detectid']

                    try:
                        j = all_detid.index(d) #match cont_tab detectid to cont.h5 detectid
                    except:
                        print(f"[{cfg.datevshot}] Warning! detectid: {d} not found in *_line.h5 file?")
                        continue

                    #No. If already marked as a continuum source, use Diagnose and need the real gmag
                    # if all_approx_gmag[j] < 640.0: #too faint (fainter than about 23 in g)
                    #     #640.0 ~ sum along the flux: flux for g = 23 @ 4640AA is about 0.32e-17 erg/s/cm2/AA
                    #     # for 0.32 * 2AA * 1000 bins = 640.0
                    #     continue

                    # rows = cont_h5.root.Spectra.read_where("detectid==d")
                    # if len(rows) != 1:
                    #     print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                    #     continue
                    # row=rows[0]
                    # gmag, gmag_unc, *_ = SU.get_best_gmag(row['spec1d']*1e-17,row['spec1d_err']*1e-17,G.CALFIB_WAVEGRID)
                    # cont_tab['gmag'][i] = gmag

                    # spec1d = cont_h5.root.Spectra.read_where("detectid==d",field="spec1d")
                    # if len(spec1d) != 1:
                    #     print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                    #     continue
                    #
                    # spec1d_err = cont_h5.root.Spectra.read_where("detectid==d",field="spec1d_err")
                    #
                    # gmag, gmag_unc, *_ = SU.get_best_gmag(spec1d[0]*1e-17,spec1d_err[0]*1e-17,G.CALFIB_WAVEGRID)

                    gmag, gmag_unc, *_ = SU.get_best_gmag(all_spec1d[j] * 1e-17, all_spec1d_err[j] * 1e-17,
                                                          G.CALFIB_WAVEGRID)
                    cont_tab['gmag'][i] = gmag

                #update
                cont_tab.write(name, format="ascii", overwrite=True)
                print(f"[{cfg.datevshot}] Updated lines table with gmag: {os.getcwd()}/{name}")

            totN = len(cont_tab)
            DETECTID = np.array(cont_tab['detectid'])
            RA = np.array(cont_tab['ra'])
            DEC = np.array(cont_tab['dec'])
            GMAG = np.array(cont_tab['gmag'])


            RMAG = np.zeros((totN,))
            IMAG = np.zeros((totN,))
            ZMAG = np.zeros((totN,))
            YMAG = np.zeros((totN,))
            SN = np.zeros((totN,))
            NAME = np.zeros((totN,))  # np.array(source_spectra['shotid'][sel])

            # re-collect the spec for those selected
            spec_2D = []
            error_2D = []
            apcor = []

            # for i in range(len(cont_tab)):
            #     d = cont_tab[i]['detectid']
            #
            #     rows = cont_h5.root.Spectra.read_where("detectid==d")
            #     if len(rows) != 1:
            #         print(f"[{cfg.datevshot}] Diagnose preselection spectra failure for {d}")
            #         spec_2D.append(np.zeros(1036))
            #         error_2D.append(np.zeros(1036))
            #         apcor.append(np.zeros(1036))
            #         continue
            #
            #     row = rows[0]
            #
            #     spec_2D.append(row['spec1d'])
            #     error_2D.append(row['spec1d_err'])
            #     apcor.append(row['apcor'])

            for i in range(len(cont_tab)):
                d = cont_tab[i]['detectid']
                try:
                    j = all_detid.index(d)  # match cont_tab detectid to cont.h5 detectid
                except:
                    #already checked in previous loop, so this should neverl trigger
                    print(f"[{cfg.datevshot}] Warning! detectid: {d} not found in *_cont.h5 file?")
                    continue

                spec_2D.append(all_spec1d[j] * dust_corr)
                error_2D.append(all_spec1d_err[j] * dust_corr)
                apcor.append(all_apcor[j])

            SPEC = np.array(spec_2D)
            ERROR = np.array(error_2D)
            WEIGHT = np.array(apcor)

            cont_h5.close()

            num_spectra = len(SPEC)

            # now write out the fits file needed for Diagnose
            outname = 'diagnose_spectra_cont.fits'
            T = Table([DETECTID, RA, DEC, NAME, GMAG, RMAG, IMAG, ZMAG, YMAG, SN],
                      names=['detectid', 'RA', 'Dec', 'shotid', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 'sn'])
            fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(T), fits.ImageHDU(SPEC),
                          fits.ImageHDU(ERROR), fits.ImageHDU(WEIGHT)]).writeto(outname, overwrite=True)

            # clean up
            del T
            del DETECTID
            del RA
            del DEC
            del GMAG
            del WEIGHT
            del RMAG
            del IMAG
            del ZMAG
            del YMAG
            del SN
            del NAME
            del SPEC
            del ERROR

            # call Diagnose' needs sklearn

            # Diagnose path is: {cfg.scriptdir}/Diagnose/diagnose.py
            print(f"[{cfg.datevshot}] Running Diagnose on {num_spectra} continuum sources ...")
            cmd = f"python3 {cfg.scriptdir}/Diagnose/diagnose.py"
            cmd += f" {outname} -li 0 -hi {num_spectra} -s 1 -q"

            system_command(cfg, cmd)

    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Exception in Diagnose.",traceback.format_exc())

    return rc


def diagnose_output_to_table(cfg):
    """

    mostly lifted from Erin's notebook ("add_z_diagnose")

    :param cfg:
    :return:
    """


    try:

        print(f"[{cfg.datevshot}] Converting Diagnose output to table...")
        os.chdir(os.path.join(cfg.cwd))

        rc = 0
        i = 0

        colnames = ['detectid', 'RA', 'Dec', 'shotid', 'gmag', 'chi2_star',
                    'chi2_galaxy', 'chi2_qso', 'z_star', 'z_galaxy', 'z_qso', 'z_best',
                    'classification', 'stellartype']
        format_dict = {'RA': '%0.6f', 'Dec': '%0.6f', 'gmag': '%0.3f', 'z_star': '%0.8f', 'z_galaxy': '%0.5f',
                       'z_qso': '%0.5f', 'chi2_star': '%0.3f', 'chi2_galaxy': '%0.3f',
                       'chi2_qso': '%0.3f', 'z_best': '%0.8f'}

        #_000.fits are the emission line dets, _001.fits are the continuum sources
        filenames = [
            "classification_000.fits",
            "classification_001.fits"
            ]

        try:
            N = fits.open(filenames[0])[2].shape[0] + fits.open(filenames[1])[2].shape[0]
            t = fits.open(filenames[0])[6].data
        except:
            print(f"[{cfg.datevshot}] Error! diagnose_output_to_table : cannot read {filenames[0]}")
            return -1

        detectid, RA, Dec, shotid, gmag, rmag = [[t[name][0]] * N for name in
                                                 ['detectid', 'RA', 'Dec', 'shotid', 'gmag', 'rmag']]

        chi2_star, chi2_galaxy, chi2_qso = ([1.] * N, [1.] * N, [1.] * N)
        z_star, z_galaxy, z_qso, z_best = ([0.] * N, [0.] * N, [0.] * N, [0.] * N)
        classification, stellartype = (['GALAXY'] * N, [''] * N)

        del t

        for fn in filenames:
            print(f'[{cfg.datevshot}] : Working on %s' % fn)
            f = fits.open(fn)
            l = f[1].shape[0]

            #print(i, i + l)
            detectid[i:i + l] = f[6].data['detectid']
            RA[i:i + l] = f[6].data['RA']
            Dec[i:i + l] = f[6].data['Dec']
            shotid[i:i + l] = f[6].data['shotid']
            gmag[i:i + l] = f[6].data['gmag']

            chi2_star[i:i + l] = f['chi2'].data[:, 0]
            chi2_galaxy[i:i + l] = f['chi2'].data[:, 1]
            chi2_qso[i:i + l] = f['chi2'].data[:, 2]
            f['zs'].data[:, 0] = f['zs'].data[:, 0] / 2.99798e8
            z_star[i:i + l] = f['zs'].data[:, 0]
            z_galaxy[i:i + l] = f['zs'].data[:, 1]
            z_qso[i:i + l] = f['zs'].data[:, 2]
            best_ind = np.array(f['class'].data, dtype=int) - 1
            zbest = np.ones((l,)) * -999.
            for j in np.arange(3):
                sel = np.where(best_ind == j)[0]
                zbest[sel] = f['zs'].data[sel, j]
            z_best[i:i + l] = zbest
            class_label = np.array(['UNKNOWN'] * l)
            for j, label in zip(np.arange(4), ['STAR', 'GALAXY', 'QSO', 'UNKNOWN']):
                sel = np.where(best_ind == j)[0]
                class_label[sel] = label
            classification[i:i + l] = class_label

            try:
                stellartype[i:i + l] = np.array(f[5].data).astype(str)
            except:
                stellartypes = []
                for item in np.array(f[5].data):
                    try:
                        stellartypes.append(item[0].astype(str))
                    except:
                        #print(item[0])
                        stellartypes.append(' ')
                stellartype[i:i + l] = stellartypes
            # print(np.unique( f[5].data))
            i += l
            # print(i, detectid[-1])

        Td = Table([detectid, RA, Dec, shotid, gmag, chi2_star,
                    chi2_galaxy, chi2_qso, z_star, z_galaxy, z_qso, z_best,
                    classification, stellartype], names=colnames)

        Td.write("diagnose_classifications.tab",format="ascii",overwrite=True)

        print(f"[{cfg.datevshot}] Wrote: diagnose_classifications.tab")

    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Exception in Diagnose_output_to_table.", traceback.format_exc())


    return rc


def prep_elixer(cfg):
    """

    prepare, but do not run, elixer for line and continuum sources

    :param cfg:
    :return:
    """


    try:
        rc = 0
        os.chdir(os.path.join(cfg.cwd))

        print(f"[{cfg.datevshot}] Preparing detctions lists for ELiXer ...")
        #make subdirs elixer/line, elixer/cont
        elixdir = os.path.join(cfg.cwd,"elixer/")

        #should purge the directory if it exists? or move to another directory? (elixer.backup?)
        if os.path.exists(elixdir):
            #see if anything has actually been done with it? is there evidence of a slurm run

            slurmfiles = glob.glob(f"{os.path.join(elixdir,'out')}/ELIXER.o*")
            elixerfiles = glob.glob(f"{os.path.join(elixdir, 'out')}/*.h5")
            if len(slurmfiles) != 0 or len(elixerfiles) != 0:
                fns = glob.glob(f"{os.path.join(cfg.cwd,'elixer')}.*")
                ct = len(fns) + 1
                outdir = os.path.join(cfg.cwd, f"elixer.backup_{ct}")
                while os.path.exists(outdir):
                    ct += 1
                    outdir = os.path.join(cfg.cwd,f"elixer.backup_{ct}")
                #now move
                print(f"[{cfg.datevshot}] Old elixer directory found. Moving it to: {outdir}")
                Path(elixdir).rename(Path(outdir))
            else: #wipe out the directory and start over
                print(f"[{cfg.datevshot}] Old elixer directory found. Unused. Will remove and rebuild.")
                shutil.rmtree(elixdir)

        Path(elixdir).mkdir(parents=True, exist_ok=True)

        if FilterDetsOnBadAmps:
            if not FORCE_CONTINUE:
                try:
                    #we are in the shot working dir
                    h5 = tables.open_file(f"{cfg.datevshot}.h5",mode='r')
                    bad_amps_list = list(h5.root.AmpStats.read_where("flag==0", field="multiframe").astype(str))
                    print(f"[{cfg.datevshot}] Loaded {len(bad_amps_list)} bad amps ...")
                    h5.close()
                except:
                    print(f"[{cfg.datevshot}] Could not load bad_amps_list")
                    bad_amps_list = []
            else:
                print(f"[{cfg.datevshot}] --force specified. Will not restrict bad amps.")
                bad_amps_list = []
        else:
            bad_amps_list = []


        make_lines = True
        tab = Table.read(f"{cfg.datevshot}_line_sourcecat.tab",format="ascii")
        line_ct = -1
        if len(tab) == 0:
            print(f"[{cfg.datevshot}] Error! no line sources recorded. ")
            make_lines = False
        else:
            sel = np.array([x['multiframe'] not in bad_amps_list for x in tab])
            print(f"[{cfg.datevshot}] Excluding {len(sel)-np.count_nonzero(sel)} / {len(sel)} line detections as residing on bad amps.")
            line_dets = list(tab['detectid'][sel])
            line_ct = len(line_dets)
            np.savetxt(os.path.join(elixdir, "line.dets"), line_dets, fmt="%d")
            print(f"[{cfg.datevshot}] Wrote {len(line_dets)} emission line detections for ELiXer to examine.")
            del tab

        make_conts = True
        tab = Table.read(f"{cfg.datevshot}_cont_sourcecat.tab",format="ascii")
        cont_ct = -1
        if len(tab) == 0:
            print(f"[{cfg.datevshot}] Error! no continuum sources recorded. ")
            make_conts = False
        else:
            sel = np.array([x['multiframe'] not in bad_amps_list for x in tab])
            print(f"[{cfg.datevshot}] Excluding {len(sel) - np.count_nonzero(sel)} / {len(sel)} continuum detections as residing on bad amps.")
            cont_dets = list(tab['detectid'][sel])
            cont_ct = len(cont_dets)
            np.savetxt(os.path.join(elixdir, "cont.dets"), cont_dets, fmt="%d")
            print(f"[{cfg.datevshot}] Wrote {len(cont_dets)} continuum detections for ELiXer to examine.")
            del tab


        tasks_per_node = 0 #use the default 40

        if line_ct > 500 or cont_ct > 500:
            tasks_per_node = 32 #slow it down and conserve memory

        if make_lines or make_conts:
            shot_h5 = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")
            line_h5 = os.path.join(cfg.cwd,f"{cfg.datevshot}_line.h5")
            cont_h5 = os.path.join(cfg.cwd, f"{cfg.datevshot}_cont.h5")
            diagnose_tab = os.path.join(cfg.cwd,f"diagnose_classifications.tab")

            which_elixer = f"python {elixer_path}/selixer.py" #"selixer.test "
            merge_name = f"elixer_{cfg.datevshot}_cat.h5"
            #get the -A ASTXXXXX from the environment and pass it to --alloc for elixer
            slurm_alloc = None
            try:
                slurm_alloc = os.getenv('SLURM_TACC_ACCOUNT')
            except:
                print(f"[{cfg.datevshot}] Could not get SLURM allocation account.")
            elixer_base_cmd = f" -f --slurm 0 --nodes 1 --log info --shot_h5 {shot_h5} --diagnose {diagnose_tab} " \
                              f" --png --error 3.0 --neighborhood 10.0 --ntasks_per_node {tasks_per_node} " \
                              f" --timex 0.8 --post_merge 2  --merge_name {merge_name} --cnn"
                            #reduce time by about 20% ... the neighborhood is faster since only checking THIS shot
            if slurm_alloc is not None:
                elixer_base_cmd += f" --alloc {slurm_alloc}"


            if len(cfg.email) > 5:
                elixer_base_cmd += f" --email {cfg.email} "

            #combine into a single job at the end, so their names need to match
            elixer_line_cmd = f" --name out --dets line.dets  --hdf5 {line_h5} "
            elixer_cont_cmd = f" --name out --continuum --dets cont.dets  --hdf5 {cont_h5} "

            minutes = 0
            dispatch_base = 0
            os.chdir(elixdir)

            #get count of line detects and continuum detects
            try:
                with open("line.dets","r") as f:
                    line_ct = np.sum([1 for line in f])
            except:
                line_ct = None

            try:
                with open("cont.dets","r") as f:
                    cont_ct = np.sum([1 for line in f])
            except:
                cont_ct = None

            # todo: act if counts are abnormally large?
            # tasks_per_node = 56 #normal default
            # if very many continuum (more than 500 or more than lines) could be an issue with memory and need to cut down
            #


            if make_lines:
                with open(os.path.join("elixer_line_cmd.txt"),"w") as f:
                    f.write(f"{which_elixer} {elixer_base_cmd} {elixer_line_cmd} \n")

                print(f"[{cfg.datevshot}] Preparing ELiXer (line) SLURM ... ")
                system_command(cfg, "source elixer_line_cmd.txt")

                try:
                    #get the time to add to the continuum sources AND get the dispatch_base
                    if safe_cd(os.path.join(elixdir,"out")):
                        fns = sorted(glob.glob("dispatch_*"))
                        if len(fns) > 0:
                            dispatch_base = int(fns[-1][-4:]) + 1 #need the +1 so you start conts at the next dispatch
                            # no point in getting the time if this failed, so assuming it was successful:
                            with open("elixer.slurm",'r') as slurm_file:
                                for line in slurm_file:
                                    if "#SBATCH -t" in line and "Run time" in line:
                                        if "day" in line: #this is a problem (could be day or days)
                                            #example: #SBATCH -t 2 days, 1:54:00            # Run time (hh:mm:ss)

                                            print(f"[{cfg.datevshot}] WARNING!!! Excessive ELiXer run time {line.rstrip()}")
                                            print(f"[{cfg.datevshot}] You may need to manually configure.")

                                            toks = line.split()

                                            t_idx = toks.index("-t")
                                            t_idx += 1 #now on the integer leading "days,"

                                            days_as_minutes = 1440 * int(toks[t_idx])

                                            t_idx += 2
                                            t = datetime.strptime(toks[t_idx], "%H:%M:%S")
                                            minutes = int(t.hour * 60 + t.minute) + days_as_minutes
                                        else:
                                            toks = line.split()

                                            t_idx = toks.index("-t")
                                            t_idx += 1
                                            t = datetime.strptime(toks[t_idx],"%H:%M:%S")
                                            minutes = int(t.hour * 60 + t.minute)
                                            break

                        #make a copy of the elixer.run as elixer_lines.run so can prepend
                        system_command(cfg,"cp elixer.run elixer_line.run")
                except:
                    print(f"[{cfg.datevshot}] Error bulding elixer.slurm file(s) ...", traceback.format_exc())

            os.chdir(elixdir)
            if make_conts:
                # is there a line slurm already set up?
                print(f"[{cfg.datevshot}] Extra time ({minutes}) and dispatch start ({dispatch_base})")
                if minutes > 0 and dispatch_base > 0:
                    elixer_cont_cmd += f" --time_add {minutes}  --dispatch_base {dispatch_base} "

                with open(os.path.join("elixer_cont_cmd.txt"), "w") as f:
                    f.write(f"{which_elixer} {elixer_base_cmd} {elixer_cont_cmd} \n")

                print(f"[{cfg.datevshot}] Preparing ELiXer (cont) SLURM ... ")
                system_command(cfg, "source elixer_cont_cmd.txt")

                # make a copy of the elixer.run as elixer_lines.run so can prepend
                if safe_cd(os.path.join(elixdir,"out")):
                    if os.path.exists("elixer_line.run"):
                        system_command(cfg, "mv elixer.run elixer_cont.run")
                        #for some reason, a simple cat is not working ....
                        #system_command(cfg, "sleep 5; cat elixer_line.run elixer_cont.run > elixer.run")
                        with open("elixer.run","w") as elixer_run:
                            with open("elixer_line.run","r") as elixer_line:
                                for line in elixer_line:
                                    elixer_run.write(line)
                            with open("elixer_cont.run", "r") as elixer_cont:
                                for line in elixer_cont:
                                    elixer_run.write(line)

            #print(f"[{cfg.datevshot}]Preparing ELiXer SLURMS ... ")
            #os.chdir(elixdir)
            #system_command(cfg,"source elixer_cmd.txt")

            print("Default ELiXer SLURMS prepared. However, you may edit ./elixer/elixer_*_cmd.txt re-run if needed.")
            print("NOTICE!!! if there is a single 'out' directory, then the line and cont files have been merged together.")
            print("To re-run, remove the ./elixer/line and ./elixer/cont directories then 'source' the edited 'elixer_cmd.txt' file."
                  " It will set up, but not queue, two SLURM jobs (line and cont).")
            print("You may then edit ./elixer/line/elixer.slurm and ./elixer/cont/elixer.slurm if needed and sbatch when ready.")
        else:
            print(f"{cfg.datevshot} prep_elixer: no sources ")

    except:
        print(traceback.format_exc())
        rc = -1
        print(f"Exception in prep_elixer():  {cfg.datevshot}")

#

########################################################################
########################################################################
########################################################################
# Main (section)
#   notice: no actual main function
########################################################################
########################################################################
########################################################################

###########
# setup
###########


if cfg.update_only:
    update_only(cfg)
    exit(0)

if queue_elixer:
    run_queue_elixer(cfg)
    exit(0)

if prep_compress > -1:
    run_prep_compress(cfg,max_simultaneous=prep_compress)
    exit(0)

if cfg.clean_only:
    print(f"[{cfg.datevshot}] Performing only the CLEAN, level : {cfg.clean} ...")
    post_clean(cfg)
    Quit(cfg,0,"Clean complete. Exiting",do_write_status=False,do_post_clean=False,do_write_summary=False) #just ran post_clean, don't need it twice

rc = precheck(cfg)
if rc < 0:
    Quit(cfg,rc,"FATAL! Precheck failed. Reduction cannot run.",do_write_status=False,do_post_clean=False,do_write_summary=False)

#if not cfg.multifits_only:
cfg.numexp, cfg.gettar_fn = num_exposures_in_shot(cfg.shotid)

do_initial_setup = True
if cfg.numexp <= 0: # and not cfg.multifits_only:
    # might not be fatal ... could be a newer observation that is not in the gettars yet, but is still accessible
    #Quit(cfg, -1, f"FATAL! Could not find shot {cfg.datevshot}",do_write_status=False)
    print(f"[{cfg.datevshot}] Did not find shot in standard gettar location. Not necessarily fatal. Will attempt to proceed ...")
    #try the initial setup anyway
    rc = initial_setup(cfg)
    do_initial_setup = False
    if rc < 0:
        Quit(cfg, rc, "Could not complete initial setup.", do_write_status=False, do_write_summary=False)
    get_local_ra_dec(cfg)  # base info, but most will be overwritten as we go
    #now we need to build the local gettar
    cfg.numexp, cfg.gettar_fn = make_local_gettar_file(cfg)
    if cfg.gettar_fn is None:
        #now this is fatal
        Quit(cfg, -1, f"FATAL! Could not find shot {cfg.datevshot}", do_write_status=False, do_write_summary=False)



# if cfg.simul == 1:
#     NumProcs_mp_rcal = 10
#     NumProcs_mp_rf1 = 10

if cfg.exp <= 0 and cfg.multifits_only:
    pass #this is fine
elif cfg.exp <= 0 and not cfg.multifits_only:
    print(f"[{cfg.datevshot}] Working on {cfg.datevshot} with {cfg.numexp} exposure(s) ...")
    if cfg.numexp == 1 or cfg.numexp == 3: #okay
        pass
    elif not cfg.multifits_only:
        #print(f"[{cfg.datevshot}] !!! bad !!! Unusual number of exposures ({cfg.numexp}).")
        print(f"********************************************************************************************")
        print(f"!!! WARNING !!! Unusual number of exposures ({cfg.numexp}) for {cfg.datevshot} !!! Reduction may be problematic.")
        print(f"                You may want to consider reducing each exposure individually. ")
        print(f"********************************************************************************************")
else:
    if cfg.exp <= cfg.numexp:
        print(f"[{cfg.datevshot}] Working on {cfg.datevshot} exposure #{cfg.exp} ...")
        if int(cfg.datevshot[0:8]) < 20240800 and cfg.numexp == 3:
            try:
                dex_dvs = np.loadtxt("/corral-repl/utexas/Hobby-Eberly-Telesco/detect/fwhm.all",dtype=str,usecols=0)
                if cfg.datevshot in dex_dvs:
                    print(f"WARNING! {cfg.datevshot} is a original HETDEX observation.")
                    print(f"         Attempting to use only one exposure may not generate the expected results.")
                #else: we are not going to print anything ... while this occurred during the original HETDEX run
                #      this is not a previously reduced HETDEX observation
            except: #just assume it could be an original observation and warn
                print(f"WARNING! {cfg.datevshot} appears to be an original HETDEX observation.")
                print(f"         Attempting to use only one exposure may not generate the expected results.")
    else:
        if cfg.multifits_only: #asssume the user knows what is going on and just let the requested exposure number stand
            if cfg.numexp == 0:
                cfg.numexp = cfg.exp
        else:
            Quit(cfg, -1, f"Invalid exposure. Requesting exp #{cfg.exp} but {cfg.datevshot} has only {cfg.numexp}",
                 do_write_status=False,do_write_summary=False)

if do_initial_setup:
    rc = initial_setup(cfg)

    if rc < 0:
        Quit(cfg,rc,"Could not complete initial setup.",do_write_status=False,do_write_summary=False)

    get_local_ra_dec(cfg) #base info, but most will be overwritten as we go

# *** DEBUG HERE ***
#get_multifits_file_status(cfg)

#print("Temporary ... make summary.txt ONLY")
#Quit(cfg,rc=0,msg="Temporary ... make summary.txt ONLY",do_write_status=True)


rc = node_setup(cfg)

if cfg.special == 1:
    print("*** SPECIAL (1) *** forcing guider fwhm then terminating")
    cfg.guider_fwhm = get_guider_fwhm(cfg)
    Quit(cfg, 0, f"Done with special handling. {cfg.datevshot}",do_write_status=False,do_write_summary=False)


if not cfg.multifits_only and (cfg.numexp < 3 or not cfg.hetdex_original):
    #print(f"Fewer than 3 exposures (assume dithers). Checking guider for seeing FWHM...")
    #print(f"[{cfg.datevshot}] Checking guider for seeing FWHM (this can take a while) ...")
    cfg.guider_fwhm = get_guider_fwhm(cfg) #this also gets the exposure times
    if cfg.guider_fwhm is not None:
        print(f"[{cfg.datevshot}] Using guider FWHM = {cfg.guider_fwhm}")
    else:
        print(f"[{cfg.datevshot}] Unable to obtain guider seeing FWHM. Will measure as best can be from available data.")

if not cfg.multifits_only and (cfg.total_exp_time is None or cfg.total_exp_time == 0):
    get_exposure_times(cfg)




# print(f"!!! DEBUG !!!")
# check_line_detections_by_ifu(cfg)
#
# Quit(cfg, 0, f"DEBUG exit. {cfg.datevshot}",do_write_status=False,do_write_summary=False)

#########
# after the initial setup, move stdout and stderr to a log file
#########

# print(f"**** DEBUG  collapsed image test ****")
# avg_sky = get_avg_sky(cfg)
# Quit(cfg,0,"Done Test")

#
# print(f"TESTING match_pngs... ")
#
# build_shot_h5(cfg)
#
# print(f"EXITING match_pngs")
# exit(0)



#cfg.orig_stdout = sys.stdout
#cfg.orig_stderr = sys.stderr
cfg.file_stdout = open(f"{cfg.datevshot}.log","a")
print(f"[{cfg.datevshot}] Logging redirected to: {cfg.cwd}/{cfg.file_stdout.name}")
#sys.stderr = cfg.file_stdout
#sys.stdout = cfg.file_stdout


# get the progress state. Useful if resuming (implied)
dtprog = progress_init(cfg)



#update for long exposures
if not cfg.multifits_only and not dtprog['s02_vdrp']:
    #if vdrp has already run, there is no point in updating
    update_vdrp_config_limits(cfg)

#begin
with open("status.run", "w") as f:
    f.write(f"[{cfg.datevshot}] running .... \n")

###########
# step1
###########

if s01_run1s and not dtprog["s01_run1s"]:
    run_run1s(cfg)

    # run any checks
    if check_run1s(cfg) < 0:
        Quit(cfg, -1, "FATAL. Initial extraction and/or base calibration failure.")

    #todo: this would be manual here, I think, but CAN copy /red1/xxx to /scratch/local/projects
    #  all the various CoFe*.fits and multi*.fits ... these are also in the
    #  local d<shot><exp> folder in the two tar files (_co.tar and _mu.tar for the CoFe*.fits and multi*.fits respectively)

    progress_update(cfg,dtprog,"s01_run1s")
else:
    print(f"[{cfg.datevshot}] Skipping s01_run1s run1s")

    #However, still need this later for vdrp pathing
    if not os.path.exists("./reductions"):
        system_command(cfg,"ln -s ../reductions reductions")


if s01b_amp_images and not dtprog["s01b_amp_images"]:
    if make_amp_images(cfg) < 0: #this is non-fatal if it fails
        print(f"[{cfg.datevshot}] creation of IFU+amp diagnostic images failed. Non-fatal. Will continue.")
    else:
        if make_amp_images(cfg,ratio=True) < 0:  # this is non-fatal if it fails
            print(f"[{cfg.datevshot}] creation of IFU+amp diagnostic RATIO images failed. Non-fatal. Will continue.")
        else:
            progress_update(cfg,dtprog,"s01b_amp_images")

else:
    print(f"[{cfg.datevshot}] Skipping s01b_amp_images")

if cfg.multifits_only:
    print(f"[{cfg.datevshot}] multi*fits files generated. "
          f"Check paths under ./reductions/{cfg.datevshot[0:8]}/virus/virus*{cfg.datevshot[-3]}/")
    print(f"[{cfg.datevshot}] Also look for diagnostic IFU+amp images under sci{cfg.datevshot}/analysis/ifus")
    Quit(cfg, 0, f"Done. Enforcing --multifits_only switch. {cfg.datevshot}", do_write_status=False)


#quick assignment ... may need these later
get_multifits_file_status(cfg)


###########
# step2
# VDRP
###########

if s02_vdrp and not dtprog["s02_vdrp"]:
    run_vdrp(cfg)

    if check_vdrp(cfg) < 0:
        Quit(cfg, -1, "FATAL. Could not make astrometric solution.")

    #todo: optional manual step here (need to be hetdex user), copy the *.dithall to
    #  /scratch/projects/hetdex/detect/dithall   (and /coral-repl/...)

    progress_update(cfg,dtprog, "s02_vdrp")

else:
    print(f"[{cfg.datevshot}] Skipping s02_vdrp vdrp")



###########
# step3
# detect (Flux Calibration)
# rallcal stuff
###########

if s03_fluxcal and not dtprog["s03_fluxcal"]:


    # #first need to prepare the reduction directory where the downstream code looks for the multi*fits
    if prepare_reduction_dir(cfg) < 0:
        Quit(cfg, -1, "FATAL. Could not untar multi*fits to reduction directory.")

    star_cat_list = [] #keep the order

    #this is largely for parallels and we want GAIA for that, so have it first
    # note: for HETDEX it was SDSS first for flux calibration and GAIA for astrometry

    # if cfg.hetdex:
    #     #SDSS for flux calibration and GAIA for astrometry
    #     if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.sdss")):
    #         star_cat_list.append("sdss")
    #
    #     if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.gaia")):
    #         star_cat_list.append("gaia")
    #
    #     if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.panstarrs")):
    #         star_cat_list.append("panstarrs")
    # else:
    #     if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.gaia")):
    #         star_cat_list.append("gaia")
    #
    #     if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.sdss")):
    #         star_cat_list.append("sdss")
    #
    #     if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.panstarrs")):
    #          star_cat_list.append("panstarrs")


    if cfg.starcat_cal is None or cfg.starcat_cal == 'sdss':
        if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.sdss")):
            star_cat_list.append("sdss")

        if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.gaia")):
            star_cat_list.append("gaia")
    elif cfg.starcat_cal == 'gaia':
        if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.gaia")):
            star_cat_list.append("gaia")

        if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.sdss")):
            star_cat_list.append("sdss")
    else: #shold not happen
        print(f"[{cfg.datevshot}] flux calibration: unexpected value in cfg.starcat_cal {cfg.starcat_cal}. Using SDSS default.")

        if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.sdss")):
            star_cat_list.append("sdss")

        if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.gaia")):
            star_cat_list.append("gaia")

    #regardless, panstarrs is always third
    if os.path.exists(os.path.join(cfg.cwd, f"vdrp/shifts/dithall.panstarrs")):
        star_cat_list.append("panstarrs")

    if len(star_cat_list) == 0:
        Quit(cfg, -1, "FATAL. Something wrong. No star catalogs available under vdrp/shifts.")

    print(f"[{cfg.datevshot}] attemping flux calibrations in this order: {star_cat_list} ...")
    for star_cat in star_cat_list:
        #print(f"[{cfg.datevshot}] flux calibration: {star_cat}")
        run_fluxcalibration(cfg,star_cat)

        if check_fluxcalibration(cfg) < 0:
            print(f"[{cfg.datevshot}] flux calibration: {star_cat} failed. Trying next ... ")
        else:
            break #this one was good

    #todo: optional:  copy detect/tp/yyyymmddvssssedtp_f.dat  to /scratch/projects/hetdex/detectp and /corral-repl/xxx

    #todo: optional: update  /scratch/projects/hetdex/detect/fwhm.all and norm.all
    #                see update_fwhm_norm script
    #

    #could be that none worked

    if cfg.hetdex_original:
        print(f"[{cfg.datevshot}] If replacing, update /corral and /scratch/projects   detect/tp/<datevshot>sedtp_f.dat")
        print(f"[{cfg.datevshot}] If replacing, update /scratch/projects/hetdex/detect/fwhm.all and norm.all ... see update_fwhm_norm script")

    #update the fwhm.out with the cfg.guider_fwhm, if it is set
    if cfg.guider_fwhm is not None and cfg.guider_fwhm > 0:
        try:

            fwhm_fn = os.path.join(cfg.cwd,f"detect/{cfg.datevshot}/fwhm.out")
            fwhm_cp = os.path.join(cfg.cwd,f"detect/{cfg.datevshot}/fwhm.out.original")

            #make a copy
            system_command(cfg,f"cp {fwhm_fn} {fwhm_cp}")

            out = np.loadtxt(fwhm_fn) #should just be a single line

            #have to write out as file since not all floats
            print(f"[{cfg.datevshot}] Updating measured seeing FWHM {out[0]} with guider reported FWHM {cfg.guider_fwhm}")
            with open(fwhm_fn,"w") as f:
                f.write(f"{cfg.guider_fwhm}\t{out[1]}\t{int(out[2])}\n")

        except:
            print(f"[{cfg.datevshot}] Could not update fwhm.out with guider.fwhm")
            print(traceback.format_exc())

    progress_update(cfg,dtprog,"s03_fluxcal")
else:
    print(f"[{cfg.datevshot}] Skipping s03_fluxcal flux calibration")


###########
# step4
# Sky Subtractions, flux calibrations
# getcen and more stuff
###########

if s04_make_shot and not dtprog["s04_make_shot"]:

    if s04a_get_ifucens and not dtprog["s04a_get_ifucens"]:
        print(f"[{cfg.datevshot}] Getting IFU centers ...")
        if run_make_ifucen(cfg) != 0:
            Quit(cfg,-1,"FATAL. Failed to get IFU centers.")
        else:
            progress_update(cfg,dtprog, "s04a_get_ifucens")
    else:
        print(f"[{cfg.datevshot}] Skipping s04a_get_ifucens IFU centers ...")

    if s04b_rfft and not dtprog["s04b_rfft"]:
        print(f"[{cfg.datevshot}] Running rfft (this may take a while) ...")
        if run_rfft(cfg) != 0:
            if cfg.avg_sky >= FAIL_AVG_SKY:
                Quit(cfg,-1,f"[{cfg.datevshot}] Average Sky is catastrophically large and/or unable to be fit: {cfg.avg_sky:0.1f}.")
            else:
                Quit(cfg, -1, "FATAL. rfft fail. One or more expected outputs failed.")
        else:
            progress_update(cfg,dtprog, "s04b_rfft")
    else:
        print(f"[{cfg.datevshot}] Skipping s04b_rfft rfft")


    if s04c_rcal_all and not dtprog["s04c_rcal_all"]:
        print("Running rcal_all ...")
        rc = run_rcal(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rcal_all fail.")
        elif rc > 0:
            print(f"[{cfg.datevshot}] rcal_all: Limited success. Non-fatal. Will continue")

        progress_update(cfg,dtprog, "s04c_rcal_all")
        #else keep going

    else:
        print(f"[{cfg.datevshot}] Skipping s04c_rcal_all rcal_all")

    if s04d_shot_h5 and not dtprog["s04d_shot_h5"]:
        rc = build_shot_h5(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. Could not build shot h5 file. Cannot continue with catalog creation.")
        progress_update(cfg,dtprog, "s04d_shot_h5")
    else:
        print(f"[{cfg.datevshot}] Skipping s04d_shot_h5 build of shot HDF5 file.")

    # check stats
    if s04e_amp_stats and not dtprog["s04e_amp_stats"]:
        rc = amp_stats(cfg)
        if rc < 0:
            print(f"[{cfg.datevshot}] Non-fatal. Could not compute amp stats from shot h5 file. Will continue anyway with catalog creation.")

        #need to add the index first
        rc = add_fiber_index(cfg)
        if rc < 0:
            print(f"[{cfg.datevshot}] Severe, but non-fatal. Could not add fiber index. Will continue anyway with catalog creation.")
        else:
            #then can add the mask
            rc = add_fiber_mask(cfg)
            if rc < 0:
                print(f"[{cfg.datevshot}] Non-fatal. Could not update shot h5 file with fiber level masking. Will continue anyway with catalog creation.")

        progress_update(cfg,dtprog, "s04e_amp_stats")
    else:
        print(f"[{cfg.datevshot}] Skipping s04e_amp_stats compute and appending of per-amp diagnostic info.")

    #basic shot analysis (mostly images for review)
    if s04f_analysis and not dtprog["s04f_analysis"] and not cfg.multifits_only:
        rc = shot_analyisis(cfg)
        if rc < 0:
            print(f"[{cfg.datevshot}] Non-fatal. Could not complete basic shot analysis output.")
        progress_update(cfg,dtprog, "s04f_analysis")
    else:
        print(f"[{cfg.datevshot}] Skipping s04f_analysis analysis/creation of diagnostic info.")

    if dtprog["s04a_get_ifucens"] and dtprog["s04b_rfft"] and dtprog["s04c_rcal_all"] and dtprog["s04d_shot_h5"] and dtprog["s04e_amp_stats"]:
        progress_update(cfg,dtprog, "s04_make_shot")
else:
    print(f"[{cfg.datevshot}] Skipping sky subtraction")


###########
# step5
# line detections
# continuum detections
###########


#precheck

if "-linedet_parms" not in cfg.args:
    cfg.linedet_parms = (1, 0.0)  # standard non-dithered, applies to all cases unless changed

if cfg.shot_only: #we are done
    Quit(cfg, 0, f"[{cfg.datevshot}] --shot_only specified.  Will end here.")
elif cfg.numexp == 1:
    pass
elif cfg.numexp == 3:
    rc = hetdex_dither(cfg)

    if rc == 1: #all good
        cfg.dither_configuration = 1 #standard hetdex
        if "-linedet_parms" not in cfg.args:
            cfg.linedet_parms = (3, 0.5)  # HETDEX style = (3, 0.5) #must be (<int>,<float>)
        #else it was specfied, so use it
    elif rc == 0: #not HETDEX dither (that is okay but we cannot move on to source detection)
        cfg.dither_configuration = 0  # non-standard, assume multiple exposures, but not dithered
        if FORCE_CONTINUE:
            print(f"[{cfg.datevshot}] Non-standard dither configuration, but --force flag set, so will continue.")
        else:
            post_clean(cfg)
            Quit(cfg, 0,f"[{cfg.datevshot}] ({cfg.numexp}) exposures not in HETDEX dither configuration and is not "
                 f"compatible with source detection. Will end here.")
    else: #fail case
        cfg.dither_configuration = -1
        Quit(cfg, -1,f"[{cfg.datevshot}] Fatal error checking dither configuration. Will terminate here.")
else: #not 1 and not 3
    post_clean(cfg)
    Quit(cfg,0,f"[{cfg.datevshot}] Number of exposures ({cfg.numexp}) incompatible with source detection. Will end here.")


if s05_detection and not dtprog["s05_detection"]:

    if s05b_rdet_rf1 and not dtprog["s05b_rdet_rf1"]:
        print("Running rdet_rf1 (line detection) ...")
        rc = rdet_rf1(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rdet_rf1 fail.")
        elif rc > 0:
            print("rdet_rf1: Limited success. Non-fatal. Will continue")

        progress_update(cfg,dtprog, "s05b_rdet_rf1")
    else:
        print(f"[{cfg.datevshot}] Skipping s05b_rdet_rf1 rdet_rf1 (line detection)")


    #check the line detections for excesses
    #rc = check_line_detections(cfg)
    rc = check_line_detections_by_ifu(cfg)
    if rc < 0:
        Quit(cfg, -1, "FATAL. rget_rf1 (line detections) fail.")




    if s05c_rgetmax  and not dtprog["s05c_rgetmax"]:
        print(f"[{cfg.datevshot}] Running rgetmax (continuum detection) ...")
        rc = rgetmax(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rgetmax fail.")
        elif rc > 0:
            print(f"[{cfg.datevshot}] rgetmax: Limited success. Non-fatal. Will continue")

        progress_update(cfg,dtprog, "s05c_rgetmax")
    else:
        print(f"[{cfg.datevshot}] Skipping s05c_rgetmax rgetmax (continuum detection)")

    if s05e_detection_hdf5 and not dtprog["s05e_detection_hdf5"]:
        rc = build_detection_hdf5(cfg)

    progress_update(cfg,dtprog, "s05e_detection_hdf5")

    if dtprog["s05b_rdet_rf1"] and dtprog["s05c_rgetmax"] and dtprog["s05e_detection_hdf5"]:
        progress_update(cfg,dtprog, "s05_detection")

else:
    print(f"[{cfg.datevshot}] Skipping detections")


##################################################
# step6
# combine detections
# run diagnose and elixer
# make source catalogs
#################################################

if s06_catalogs and not dtprog["s06_catalogs"]:

    print(f"[{cfg.datevshot}] Catalog creation ... ")

    if s06b_fof and not dtprog["s06b_fof"]:
        try:
            #line_tab = Table.read(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.fits"),format="fits")
            line_h5 = tables.open_file(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.h5"))
            line_tab = Table(line_h5.root.Detections.read())
            line_h5.close()

            try: #bright continuum, weak (false) line cut (linewidth is sigma so need fwhm
                original_settings = np.geterr()
                np.seterr(divide='ignore', invalid='ignore')
                norm_obs_ew = line_tab['flux'] / line_tab['continuum'] / (2.355*line_tab['linewidth'])
                norm_obs_ew[np.isnan(norm_obs_ew)] = 1.0
                norm_obs_ew[np.isinf(norm_obs_ew)] = 1.0
                np.seterr(**original_settings)
            except:
                print(f"[{cfg.datevshot}] Exception computing normed Obs EW. {traceback.format_exc()}")
                norm_obs_ew = None

            #subselect "nominal" good
            if cfg.linedet_filter == 0:
                print(f"[{cfg.datevshot}] Standard extra restriction on line detections. "
                      f"Includes SNR >= {DEFAULT_MIN_SNR_FOR_ELIXER}")
                esel = np.array(line_tab['continuum'] >= -3)
                esel = esel & np.array(line_tab['sn'] >= DEFAULT_MIN_SNR_FOR_ELIXER)
                esel = esel & np.array(line_tab['chi2'] <= 2.5)
                # this is a bit more liberal than standard HETDEX_API (1.6 and 14, I think)
                esel = esel & np.array(line_tab['linewidth'] >= 1.5) & np.array(line_tab['linewidth'] <= 16)
                esel = esel & np.array(line_tab['chi2fib'] <= 4.5)  # fairly restrictive, this is from .mc file column #19
               # esel = esel &
            else:
                esel = np.full(len(line_tab),True)
                esel = esel & np.array(line_tab['linewidth'] >= 1.5) #just the lower limit
                esel = esel & np.array(line_tab['chi2fib'] <= 4.5)
                if cfg.linedet_filter >= 1.1:
                    # note: the SNR restriction is passed to HETDEX_API in --sn_min so there
                    # should not be any in the base set that are less than that level, but just to be
                    # consistent, apply it anyway
                    esel = esel & np.array(line_tab['sn'] >= cfg.linedet_filter)
                    print(f"[{cfg.datevshot}] Limited extra restriction on line detections. "
                          f"Includes SNR >= {cfg.linedet_filter}")
                else:
                    print(f"[{cfg.datevshot}] Limited extra restriction on line detections.")

            if norm_obs_ew is not None:
                #from a bit above, try to limit the many false lines found on top of positive continuum
                try:
                    reject_obs_ew = np.array(norm_obs_ew < 0.5) * np.array(line_tab['continuum'] > 0.35)
                    esel = esel * ~reject_obs_ew
                    print(f"[{cfg.datevshot}] Removing {np.count_nonzero(reject_obs_ew)} / {len(esel)} line detections for out of range normed Obs EW.")
                except:
                    print(f"[{cfg.datevshot}] Exception applying normed Obs EW. {traceback.format_exc()}")

            line_tab = line_tab[esel]

            if line_tab is not None:
                tname = f"{cfg.datevshot}_line_sourcecat.tab"
                line_tab.write(tname, format="ascii", overwrite=True)
                print(f"[{cfg.datevshot}] Initial lines source table, {len(line_tab)} rows: {os.getcwd()}/{tname}")

                fof_3d_lines_tab = make_3d_friend_table_for_shot(line_tab, dsky_3D=6.0, dwave=4.0)
                if fof_3d_lines_tab is not None:
                    fof_2d_lines_tab = make_2d_friend_table_for_shot(line_tab, dsky_2D=3.0)

                    if fof_2d_lines_tab is not None:
                        line_cat = join(fof_2d_lines_tab, fof_3d_lines_tab, keys="detectid")
                        line_cat["wave_group_id"] = MaskedColumn(line_cat["wave_group_id"]).filled(0)
                        line_cat.rename_column("id", "source_id")

                        # add groups back into cont_tab and overwrite
                        line_tab['source_id'] = -1
                        line_tab['sel_det'] = True

                        for i in range(len(line_tab)):
                            try:
                                # find the match
                                #sel = line_cat['detectid'] == line_tab['line_detectid'][i]
                                sel = line_cat['detectid'] == line_tab['detectid'][i]
                                # should be exactly one
                                if np.count_nonzero(sel) == 1:
                                    line_tab['source_id'][i] = line_cat['source_id'][sel]
                                elif np.count_nonzero(sel) == 0:
                                    continue  # this is okay, this has no group, so stands on its own
                                else:
                                    print(
                                        f"Error! Unexpected matches ({np.count_nonzero(sel)}) for {line_tab['detectid'][i]}")
                                    continue
                            except:
                                print(traceback.format_exc())

                        # now pick the seldet
                        uniq_src = sorted(np.unique(line_tab['source_id']))
                        try:
                            uniq_src.remove(-1)  # exclude/ignore the -1
                        except:
                            pass

                        for src_id in uniq_src:
                            sel = line_tab['source_id'] == src_id
                            line_tab['sel_det'][sel] = False
                            #idx = np.argmax(line_tab['lineflux'][sel])
                            idx = np.argmax(line_tab['flux'][sel])
                            sel_detid = line_tab['detectid'][sel][idx]
                            line_tab['sel_det'][line_tab['detectid']==sel_detid]=True
                            #line_tab['sel_det'][sel][idx] = True #this IS isn column:row order, not sure why not working


                        # tname = f"{cfg.datevshot}_line.fits"
                        # line_tab.write(tname, format="fits", overwrite=True)
                        tname = f"{cfg.datevshot}_line_sourcecat.tab"
                        line_tab.write(tname, format="ascii", overwrite=True)
                        cfg.num_line_dets = len(line_tab)
                        print(f"[{cfg.datevshot}] Updated lines source table, {cfg.num_line_dets} rows: {os.getcwd()}/{tname}")

                        #esel = np.array(line_tab['sel_det'] == True)
                        #np.savetxt('elixer_line.dets',line_tab['detectid'][esel],fmt="%d")

                    else:
                        print(f"[{cfg.datevshot}] Notice! (1) Could not combine lines detections by FoF or no lines qualified to cluster.")
                else:
                    print(f"[{cfg.datevshot}] Notice! (2) Could not combine lines detections by FoF or no lines qualified to cluster.")
            else:
                print(f"[{cfg.datevshot}] Error! Could not combine lines detections. Lines table not found.")
        except:
            print(f"[{cfg.datevshot}] Error! Could not combine lines detections by FoF.")
            print(traceback.format_exc())

        try:
            #cont_tab = Table.read(os.path.join(cfg.cwd,f"{cfg.datevshot}_cont.fits"),format="fits")
            cont_h5 = tables.open_file(os.path.join(cfg.cwd, f"{cfg.datevshot}_cont.h5"))
            cont_tab = Table(cont_h5.root.Detections.read())
            if "contflux" not in cont_tab.columns:
                cont_tab["contflux"] = [np.nanmedian(cont_h5.root.Spectra.read_where("detectid==x",field="spec1d")[0][200:800])
                                        for x in cont_tab['detectid']]

            cont_h5.close()
            if cont_tab is not None:
                tname = f"{cfg.datevshot}_cont_sourcecat.tab"
                cont_tab.write(tname, format="ascii", overwrite=True)
                print(f"[{cfg.datevshot}] Initial continuum source table, {len(cont_tab)} rows: {os.getcwd()}/{tname}")

                fof_3d_cont_tab = make_3d_friend_table_for_shot(cont_tab, dsky_3D=6.0, dwave=4.0)
                if fof_3d_cont_tab is not None:
                    fof_2d_cont_tab = make_2d_friend_table_for_shot(cont_tab, dsky_2D=3.0)

                    if fof_2d_cont_tab is not None:
                        cont_cat = join(fof_2d_cont_tab, fof_3d_cont_tab, keys="detectid")
                        cont_cat["wave_group_id"] = MaskedColumn(cont_cat["wave_group_id"]).filled(0)
                        cont_cat.rename_column("id", "source_id")

                        #add groups back into cont_tab and overwrite
                        cont_tab['source_id'] = -1
                        cont_tab['sel_det'] = True

                        for i in range(len(cont_tab)):
                            try:
                                #find the match in cont_cat
                                #sel = cont_cat['detectid'] == cont_tab['cont_detectid'][i]
                                sel = cont_cat['detectid'] == cont_tab['detectid'][i]
                                #should be exactly one
                                if np.count_nonzero(sel) == 1:
                                    cont_tab['source_id'][i] = cont_cat['source_id'][sel]
                                elif np.count_nonzero(sel) == 0:
                                    continue #this is okay, this has no group, so stands on its own
                                else:
                                    #print(f"Error! Unexpected matches ({np.count_nonzero(sel)}) for {cont_tab['cont_detectid'][i]}")
                                    print(f"[{cfg.datevshot}] Error! Unexpected matches ({np.count_nonzero(sel)}) for {cont_tab['detectid'][i]}")
                                    continue
                            except:
                                print(traceback.format_exc())

                        #now pick the seldet
                        uniq_src = sorted(np.unique(cont_tab['source_id']))
                        try:
                            uniq_src.remove(-1) #exclude/ignore the -1
                        except:
                            pass

                        for src_id in uniq_src:
                            sel = cont_tab['source_id'] == src_id
                            cont_tab['sel_det'][sel] = False
                            idx = np.argmax(cont_tab['contflux'][sel])
                            sel_detid = cont_tab['detectid'][sel][idx]
                            cont_tab['sel_det'][cont_tab['detectid'] == sel_detid] = True
                            #idx = np.argmax(cont_tab['flux'][sel])
                            #cont_tab['sel_det'][sel][idx] = True

                        # tname = f"{cfg.datevshot}_cont.fits"
                        # cont_tab.write(tname, format="fits", overwrite=True)
                        tname = f"{cfg.datevshot}_cont_sourcecat.tab"
                        cont_tab.write(tname, format="ascii", overwrite=True)
                        cfg.num_cont_dets = len(cont_tab)
                        print(f"[{cfg.datevshot}] Updated continuum source table, {cfg.num_cont_dets} rows: {os.getcwd()}/{tname}")
                        #lastly sub select based on some minimum contflux ??

                        #now done elsewhere
                        #esel = np.array(cont_tab['sel_det'] == True)
                        #np.savetxt('elixer_cont.dets', cont_tab['detectid'][esel],fmt="%d")


                    else:
                        print(f"[{cfg.datevshot}] Notice! (1) Could not combine continuum detections by FoF or no continuum sources qualified to cluster.")
                else:
                    print(f"[{cfg.datevshot}] Notice! (2) Could not combine continuum detections by FoF or no continuum sources qualified to cluster.")
            else:
                print(f"[{cfg.datevshot}] Error! Could not combine continuum detections. Continuum table not found.")
        except:
            print(f"[{cfg.datevshot}] Error! Could not combine continuum detections by FoF.")
            print(traceback.format_exc())

        progress_update(cfg,dtprog, "s06b_fof")
    else:
        print(f"[{cfg.datevshot}] Skipping s06b_fof Friends of Friends clustering.")
    #end if s06b_fof

    #we now have clustered lines and continuum
    #the additional clustering performed in hetdex_api is overkill or unnecessary here, I think, as this is just
    # a single shot where hetdex_api is clustering over all the shots
    #It might even not really be necessary for the 3D and 2D clustering calls ... I think the 3D is sufficient, but
    #  for now, leave them both in and combine

    #all we REALLY care about at this point is the detectids with matching source_ids
    #we want to roll those in with the original tables and then, for each source_id just
    #  use the one with the highest SNR or lineflux? for lines and the highest continuum for cont sources
    if s06c_diagnose and not dtprog["s06c_diagnose"]:

        rc = diagnose(cfg)

        if rc != 0:
            print(f"[{cfg.datevshot}] Diagnose: Limited success. Non-fatal. Will continue")

        rc = diagnose_output_to_table(cfg)
        if rc != 0:
            print(f"[{cfg.datevshot}] Diagnose conversion: Limited success. Non-fatal. Will continue")

        progress_update(cfg,dtprog, "s06c_diagnose")
    else:
        print(f"[{cfg.datevshot}] Skipping s06c_diagnose Diagnose.")

    if s06d_elixer and not dtprog["s06d_elixer"]:

        rc = prep_elixer(cfg)

        progress_update(cfg,dtprog, "s06d_elixer")
    else:
        print(f"[{cfg.datevshot}] Skipping s06d_elixer ELiXer prep.")

    #s06_catalogs = s06_catalogs | s06b_fof | s06c_diagnose | s06d_elixer | s06e_source_cat
    if dtprog["s06b_fof"] and dtprog["s06c_diagnose"] and dtprog["s06d_elixer"]:
        progress_update(cfg,dtprog, "s06_catalogs")




##########
# DONE
##########

post_clean(cfg)

Quit(cfg, 0, f"Complete: {cfg.datevshot}")

##############################
# Manual Step Reminders
##############################

print(f"\nOPTIONAL manual INSPECTION steps. You may want to: \n"
      f"* check the log files (under each subdir) for issues.\n"
      f"* examine the focal plane images (d<date>s<shot>exp<??>/*.png)\n"
      f"* examine the dither images (if dithered)\n")

print(f"\nOPTIONAL manual COPY steps. You may want to: \n"
      f"* copy {cfg.cwd_orig}/reductions/{cfg.datevshot[0:6]}/...  to another reductions directory\n"
      f"* copy {cfg.cwd}/vdrp/shifts/dithall.gaia/*.dithall to /scratch/projects/hetdex/detect/dithall/  and /corral-repl \n"
      f"        note: you would need to be the hetdex user\n"
      f"")
