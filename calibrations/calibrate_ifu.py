"""

Update the Calibration a single IFU for a given time period (typically over one calendar month)

!!! DOES NOT INCLUDE Amp to Amp (ata) or WAVELENGTH solutions, which require the full array of IFUs to calibrate !!!!
!!! Therefore, this is ONLY to update fiber to fiber, fiber trace, etc that can be done on a single IFU !!!

This is based on the HETDEX Calibration proceedure (see notes_calibration.odt) but operates on a single IFU
  and is self-contained (e.g. this script runs the full calibration with automated checks along the way)


MEMORY -- 4 simultaneous is okay on vm-small
       -- up to 32 may be okay on normal

RUNTIME -- about 1.25 hours (70 - 75 minutes) if tar files already untarred
        -- about 2 - 2.5 hours if need to copy tar files

This follows the same basic layout as ../single_shot_reducetion/reduce_shot.py
"""


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

try:
    from filelock import FileLock
except:
    print("You need to install filelock (e.g.: pip install --user filelock) ")
    exit(-1)

import tables
from astropy.table import Table, join, Column, MaskedColumn # unique, vstack, hstack
from astropy.io import fits
# noinspection PyUnresolvedReferences
from h5tools import amp_stats as AmpStats
# noinspection PyUnresolvedReferences
import hetdex_tools.fof_kdtree as fof

# noinspection PyUnresolvedReferences
from elixer import global_config as G
# noinspection PyUnresolvedReferences
from elixer import spectrum_utilities as SU
# noinspection PyUnresolvedReferences
from elixer import utilities as Utils


#just want the path for hetdex_api (see later)
import importlib.util

import traceback
import psutil
import time
import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
plt.style.use('default')

########################################################################
# CONFIGURATION
########################################################################
EchoCmds = True  #if True echo system commands to the log
FilterDetsOnBadAmps = True # if True, do NOT pass detections that are on reported bad amps to elixer for processing
DefaultClean = 1 #clean 0 does nothing, 1 cleans script files and temporary stuff, 2,3,4,5 are increasingly agressive
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot
LastKnownFplane = "fp20240731"
ElixerSnrThresh = 4.5 #do not run elixer on line sources where the S/N < 4.5
GuiderFWHM_ALL = True #if True and using the GUIDER FWHM, use all within the observation timeframe
                      #if False, just use the two nearest in time to the end of the observation

#if we know the number of active shots, use these limits / # of active shots
MaxTotalProcs_mp_rcal = 60 #flux calibration
MaxTotalProcs_mp_rf1 = 40 #line detections, 4 is about right (on averaged) to avoid memory issues on lonestar6

#10 seems to be the max for a vm-small
MaxPerShotProcs_mp_rcal = 10 #even if more Total processes are available, limit to this for any given shot
MaxPerShotProcs_mp_rf1 = 10  #even if more Total processes are available, limit to this for any given shot

#with no information, use these limite
NumProcs_mp_rcal = 6 #flux calibration
NumProcs_mp_rf1 = 4 #line detections, 4 is about right (on averaged) to avoid memory issues on lonestar6


ScriptRepo = "/work/03261/polonius/hetdex/calibrations/calscripts.single_ifu"
LocalScriptRepo = "./local_script_repo" #useful if running multiple single shots ... can copy remotely once
                                        #then copy locally from here for each shot
                                        #set to None if you do NOT want to use a local script dir cache
                                        #  and force a copy from the main repo each time

Lock_mutex_fn = "lsr.lock"  #all instaces under this directory share this
Node_basedir = "/tmp/hetcal"
Lock_tmp_mutex_fn = "/tmp/tmp_hetcal.mutex" #all instaces on one node (common /tmp) share this
WorkDirRoot = "./"
#user specific
red1path = None # if None, will use the local (cwd) as the basepath, otherwise user can edit and specify one here
                # "/scratch/03261/polonius/red1/reductions/"

HETRaw_archive = "/corral/utexas/Hobby-Eberly-Telesco/het_raw/"
#HETDEXSurvey = "/corral/utexas/Hobby-Eberly-Telesco/hdr5/survey/survey_hdr5.h5"
HETDEXSurvey = "/scratch/projects/hetdex/hdr5/survey/survey_hdr5.h5"
HET_by_date = "/work/03946/hetdex/maverick/"
karlgettar = "/work/00115/gebhardt/maverick/gettar/"
karlfplane = "/work/00115/gebhardt/maverick/fplane/"
karlhome = "/home1/00115/gebhardt"
#tarlist_home = "/work/00115/gebhardt/lib_calib/tarlists/"
tarlist_home = karlgettar #"/work/00115/gebhardt/gettar/tarlists/"
hetdex_projects_path = "/scratch/projects/hetdex/"

hetdex_api_path = os.path.dirname(importlib.util.find_spec("hetdex_api").origin)
#there is an extra "hetdex_api" at the end that points into the lower level directory for that.
#h5tools is actually a sibling
hetdex_api_path = "/".join(hetdex_api_path.split("/")[0:-1])


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
except:
    ApproxBaseRAM = -1
    MaxSafeActiveShots = 0


########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:

    clean: int = 0  # post run clean level; 0 = do not clean
    clean_only: bool = False  # if True, do nothing except for -clean
    clean_done: bool = False  # set to True if the post_clean() has been performed
    # simul : int = 0 #number of simultaneous shots being run (e.g. tasks per node), 0 is unset
    update_local_repo: bool = False
    update_only: bool = False
    overwrite: bool = False

    resume: bool = False

    email: str = ""

    cwd_orig: str = os.getcwd()
    cwd: str = os.getcwd()
    virus_tar_path: str = None
    # red1dir: str = f"{os.getcwd()}" #e.g. <path>/red1/reductions
    # BUT "reductions" is added later, so stop at the "red1" equivalent; so is same as just cwd_orig here
    scriptdir: str = ""
    gettar_fn: str = ""  # the runs* or runt* file from karlgettar folder with the date, shot, exp data

    orig_stdout = None
    orig_stderr = None
    file_stdout = None


    special = 0  # do some special, direct edit code stuff
    ifuslot = None
    ifuid = None
    specid = None
    ifu_fqid = None
    log_id = None

    yyyymm_fnbase = None #will be the actual YYYYMM for < 202408 and other values after (202488, 202500, 202600, etc)
    rtx_filelist = []
    runx_filelist = []
    shotnum = [] #list of strings of 3 digit shot numbers
    day = [] #list of days (YYYYMMDD) that match up with the shotnum


########################################################################
# Basic user input
########################################################################

args = list(sys.argv)  # python3 map is no longer a list, so need to cast here
del args[0]  # args.pop(0) #remove THIS file
args = [x.replace("--", "-") for x in args]

cfg = Config()

if "-help" in args:
    # do NOTHING else except print the help
    help = """
    
    !!! NOTICE !!!
    !!! This DOES NOT INCLUDE Amp to Amp (ata) or WAVELENGTH solutions, which require the full array of IFUs to calibrate !!!
    !!! Therefore, this is ONLY to update fiber to fiber, fiber trace, etc that can be done on a single IFU !!!

    usage: python reduce_shot.py <shot> [switches]
           where <shot> is YYYYMMDDSSS or YYYYMMDDvSSS

    switches (all optional):
    --clean <integer> : todo clean help

    --clear_mutex : Special operation - clears the concurrency mutex, resetting allowed concurrent active processes.
                    This takes priority over all other switches and will terminate with this action.

    --email <str> : if provided will attach to the elixer slurm job so this email address will get notifications

    --help : display this help text and exit

    # --ifuslot <str(3)> : work with this IFUSlot. If there is only one matching IFU for the selected date range
    #                      this is sufficent
                         
                         
    --ifustr <str(11)> : work with this fully qualified IFU string as: <specid>_<ifuslot>_<ifuid>
                         note: the normal shorthand is to refer to the ifuslot (2nd value)
                         
    --nolimit : if present, overrides the in-code limiting of simultaneous active shots per node

    --overwrite : removes the shot working directory completely and (re)starts fresh. 
              !!! Notice: --resume has priority over --overwrite


    --resume : (re)starts roughly at the last completed step (see sciXXXX/progress.dat)
           !!! Notice: This does NOT re-run steps that completed with failures, it only re-runs incomplete steps.
           !!! Notice: --resume has priority over --overwrite


    --update : removes and re-fetches the local_script_depo prior to running
               on a --resume, also updates the scripts already in the shot working directory
               
    --yyyymm : the year and month to calibrate, 6 digits

    """

    print(help)
    exit(0)

len_args = len(args)
queue_elixer = False
prep_compress = 0  # not just a boolean, use as the max simultaneous shots to process

if "-clear_mutex" in args:
    print("Hidden switch: clearing mutex, resetting allowed concurrent active processes.")
    print("This takes priority over all other switches and will terminate with this action.")
    fns = glob.glob(f"{Node_basedir}/*.sync")
    print(f"Clearing {len(fns)} active sync files.\n{[os.path.basename(x) for x in fns]}")
    os.system(f"rm {Node_basedir}/*.sync")
    print("Done. Exiting.")
    exit(0)


if "-clean" in args:
    i = args.index("-clean")
    try:
        cfg.clean = int(args[i + 1])
        if cfg.clean < 0:  # negative values are the same, but force just the clean operation (e.g. to be used on an old run)
            cfg.clean *= -1
            cfg.clean_only = True
    except:
        print(f"Invalid -clean specified")
        exit(-1)

    del args[i + 1]  # args.pop(0) #remove THIS file
    args.remove("-clean")
else:
    cfg.clean = DefaultClean  # usually this is level 1



if "-update" in args:
    cfg.update_local_repo = True
    args.remove("-update")

if "-overwrite" in args:
    cfg.overwrite = True
    args.remove("-overwrite")


if "-resume" in args:  # opposide of --overwite ... do NOT touch the (intermediate) output of the working directory
    cfg.resume = True
    args.remove("-resume")


if "-email" in args:
    i = args.index("-email")
    try:
        cfg.email = args[i + 1]
        # really basic sanity check
        if '@' not in cfg.email:
            print(f"Invalid -email format specified: {cfg.email}")
            exit(-1)
    except:
        print(f"Invalid -email specified: {args[i + 1]}")
        exit(-1)

    del args[i + 1]  # args.pop(0) #remove THIS file
    args.remove("-email")


if "-nolimit" in args:
    MaxSafeActiveShots = 0
    print("*** --nolimit Override. Do NOT enforce simultaneous active shot limit per node.")
    args.remove("-nolimit")

if "-special" in args:
    i = args.index("-special")
    cfg.special = int(args[i + 1])
    del args[i + 1]
    args.remove("-special")
    print(f"*** --special condition invoked. Condition = {cfg.special}")

# this would usually be enough to figure out the whole IFU string, but let's keep it explicit
#
# if "-ifuslot" in args:
#     i = args.index("-ifuslot")
#     try:
#         cfg.ifuslot = args[i + 1]
#         if len(cfg.ifuslot) != 3 and 1 <= int(cfg.ifuslot) <= 999:
#             print(f"Invalid -ifuslot specified: {args[i + 1]}")
#             exit(-1)
#     except:
#         print(f"Invalid -ifuslot specified: {args[i + 1]}")
#         exit(-1)
#
#     del args[i + 1]  # args.pop(0) #remove THIS file
#     args.remove("-ifuslot")

if "-ifustr" in args:
    i = args.index("-ifustr")
    try:
        ifustr = args[i + 1]
        if len(ifustr) != 11:
            okay = False
            #sanity check, could be a copy past with "multi_" prefixed
            if 17 <= len(ifustr) <= 20:
                toks = ifustr.split("_")
                if toks[0].lower() == "multi":
                    if len(toks) >= 5: #has an amp at the end if 5 or not if 4
                        ifustr = "_".join(toks[1:4])
                        okay = True

            if not okay:
                print(f"Invalid -ifustr specified: {args[i + 1]}")
                exit(-1)

        #spec slot ifuid
        toks = ifustr.split("_")
        if len(toks) != 3:
            print(f"Invalid -ifustr specified: {args[i + 1]}")
            exit(-1)

        cfg.specid = str(int(toks[0])).zfill(3)
        cfg.ifuslot = str(int(toks[1])).zfill(3)
        cfg.ifuid = str(int(toks[2])).zfill(3)
        cfg.ifu_fqid = f"{cfg.specid}_{cfg.ifuslot}_{cfg.ifuid}" #should be the same as the input string

    except:
        print(f"Invalid -ifustr specified: {args[i + 1]}")
        exit(-1)

    del args[i + 1]  # args.pop(0) #remove THIS file
    args.remove("-ifustr")

if "-yyyymm" in args:
    i = args.index("-yyyymm")
    try:
        cfg.yyyymm = args[i + 1]
        if len(cfg.yyyymm) != 6 and not (2017 <= int(cfg.yyyymm[0:4]) <= 2050) and not (1 <= int(cfg.yyyymm[-2:]) <= 12):
            print(f"Invalid -yyyymm specified: {args[i + 1]}")
            exit(-1)
    except:
        print(f"Invalid -yyyymm specified: {args[i + 1]}")
        exit(-1)

    del args[i + 1]  # args.pop(0) #remove THIS file
    args.remove("-yyyymm")

# whatever is left should be the shot
if len(args) > 0:
    #extra, unexpected arguments??
    print(f"Fatal: Unexpected additional args: {args}")
    print(f"exititing....")
    exit(-1)


#set a simple log identifier
cfg.log_id = f"ifu{cfg.ifuslot}_{cfg.yyyymm}"

print(f"[{cfg.log_id}] Initializing calibration ...")


########################################################################
# worker functions
########################################################################

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
        if os.getcwd()[-1 * len(path):] == path:
            return True
        else:
            os.chdir(orig)
            print(f"!!!!! WARNING. Failed to cd to {path} ")
            return False
    except:
        return False



def precheck(cfg):
    """
    sanity check files, directories to see if this reduction can run
    (e.g. if the lib_calib path is not available for the date, then this cannot run)

    :param cfg:
    :return:
    """

    try:

        rc = 0
        # echo a few key paths:
        print(f"[{cfg.log_id}] Precheck. HETDEX_API path: {hetdex_api_path}")


        month = cfg.yyyymm[0:6]
        path_check = os.path.join(hetdex_projects_path,f"lib_calib/{month}")
        if not os.path.exists(path_check):
            #todo: this might not be fatal ... could calibrate anyway??  Should think about that.
            # would we automatically copy the final output here? Needs to be a manual step, since
            # must be the hetdex user to do the copy
            print(f"Precheck fail. HETDEX lib_calib path does not exist: {path_check}")
            return -1

        if not os.path.exists(ScriptRepo):
            print(f"Precheck fail. ScriptRepo does not exist: {ScriptRepo}")
            return -1

        #common missing installs (that don't show up until later)


        if rc != 0:
            return -1

    except:
        print(f"Exception in precheck: {traceback.format_exc()}")

    return 0


def get_calibration_tarlist(cfg):
    """
    get the tar list for calibration
    depends on the ifuslot and the yyyymm

    e.g. rt1.202600 rt1b.202600  rt1c.202600 rt2.202600 rt3.202600 rta.202600
    e.g and runt202600 and runs202600

    :param cfg:
    :return:
    """

    rt = []
    run = []
    fnbase = None
    try:
        if cfg.yyyymm <= "202407": #201701 to 202407
            fnbase = f"{cfg.yyyymm}"
        elif cfg.yyyymm <= "202501":
            fnbase = "202488" #should be 2024xx??
        else: #2025, 2026, etc
            fnbase = f"{cfg.yyyymm[0:4]}00"

        rt = sorted(glob.glob(os.path.join(karlgettar,f"rt*.{fnbase}")))
        run = sorted(glob.glob(os.path.join(karlgettar, f"run?{fnbase}")))  #runt and runs

    except:
        print(f"Exception in get_calibration_tarlist(): {traceback.format_exc()}")

    return fnbase, rt, run


def initial_setup(cfg):
    """
    copy from script repo(s)

    change the cwd to the new workdir for this shot

    this is largely equivalent to what rsetups would do, but here for a single shot and not a whole month

    :return:
    """

    #check that the tarlist exists
    # this is more complicated and is done later
    # tarlist_fqfn = os.path.join(tarlist_home,f"{cfg.yyyymm}tarlist")
    # if not os.path.exists(tarlist_fqfn):
    #     print(f"FATAL. Could not find tarlist: {tarlist_fqfn}")
    #     return -1

    #setup temporary het_raw storage
    if not os.path.exists(f"{cfg.cwd_orig}/het_raw"):
        os.makedirs(f"{cfg.cwd_orig}/het_raw", exist_ok=True)

    #setup local outputdir
    cfg.lib_calib = os.path.join(os.getcwd(),f"lib_calib/{cfg.yyyymm}")
    os.makedirs(cfg.lib_calib,exist_ok=True)
#>>>> HERE <<<<

    workdir = os.path.join(WorkDirRoot,f"cal{cfg.yyyymm}/multi_{cfg.ifu_fqid}")

    resume = False #notice: cfg.resume MAY be true and this can still be false if the directory does not already exist
    if os.path.exists(workdir):
        if cfg.resume:
            print(f"[{cfg.log_id}] Resuming. Leave directory intact: {workdir}")
            resume = True
        elif cfg.overwrite:
            print(f"[{cfg.log_id}] Overwriting directory {workdir} ... ")
            shutil.rmtree(workdir)
        else:
            print(f"[{cfg.log_id}] Calibration directory already exists here! {workdir}")
            print(f"[{cfg.log_id}] Please include --resume or --overwrite to make intention clear.")
            print(f"[{cfg.log_id}] Or is this a repeated SLURM task?")
            if cfg.clean > 1:
                print(f"[{cfg.log_id}] Did you intend to only clean up the directory? "
                      f"If so, --clean needs to be the negative value. (see --help)")
            return -1

    if not resume:
        os.makedirs(workdir)

        if LocalScriptRepo is not None:
            #need to see if can get the file lock in case another instance is copying
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not calXXXX)
            with lock:
                if cfg.update_local_repo:
                    print("Updating local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, ""),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print("Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    #copy first to local script repo
                    print("Copying to local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, ""),
                                    os.path.join(os.getcwd(),LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

            #lock auto releases
        else:
            print("Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"")
    else:
        if LocalScriptRepo is not None:
            fatal_rtn = False
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:

                if cfg.update_local_repo:
                    print("Updating local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, ""),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print("Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    print("Fatal! --resume selected, but no script repo.")
                    fatal_rtn = True
            # lock auto releases
            if fatal_rtn:
                return -1
        else:
            print("Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"")


    os.chdir(workdir)
    cfg.cwd = os.getcwd() #now under the  cal<yyyymm>/ifu_xxx directory


    if not resume or cfg.update_local_repo:

        # some other process might be updating the local repo ... do not proceed IF there is a lock
        lock = FileLock(Lock_mutex_fn)  # we are in the top directory (not sciXXXX)
        with lock:
            # obtained the lock, so we should be good now,
            # just release and go
            print(f"[{cfg.log_id}] Mutex checked. Okay. Safe to local_repo.")
            # lock auto releases

        print(f"Copying source code to working directory {cfg.cwd}...")

        ## if ANY of this fails it is fatal
        #for rstep1
        shutil.copy2(os.path.join(cfg.scriptdir, "run1t"), ".") #needs path edits, cp call edits
        system_command(cfg, f"sed -i s#runtChangeMe#runt{cfg.yyyymm}# run1t")
        system_command(cfg, f"sed -i s#scriptdir_ChangeMe#{cfg.cwd}/# run1t")
        system_command(cfg, f"sed -i s#het_raw_ChangeMe#{cfg.cwd_orig}/het_raw# run1t")



        shutil.copy2(os.path.join(cfg.scriptdir, "rback"), ".") #needs path edits, cp call edits
        system_command(cfg, f"sed -i s#scriptdir_ChangeMe#{cfg.cwd}# rback")
        #rj is built during rstep1
        shutil.copy2(os.path.join(cfg.scriptdir, "rbfits"), ".")
        system_command(cfg, f"sed -i s#het_raw_ChangeMe#{cfg.cwd_orig}/het_raw# rbfits")
        shutil.copy2(os.path.join(cfg.scriptdir, "rbfits0"), ".")
        system_command(cfg, f"sed -i s#het_raw_ChangeMe#{cfg.cwd_orig}/het_raw# rbfits0")
        shutil.copy2(os.path.join(cfg.scriptdir, "rimarb"), ".") #may be fine as is
        shutil.copy2(os.path.join(cfg.scriptdir, "rback1"), ".") #may need path edits , cp call edits ?
        shutil.copy2(os.path.join(cfg.scriptdir, "sun_use.dat"), ".") #fine as is


        #for rstep2
        shutil.copy2(os.path.join(cfg.scriptdir, "rgetcal0"), ".") #need to edit cp calls
        shutil.copy2(os.path.join(cfg.scriptdir, "rgetcmb1"), ".") #may be okay as is
        shutil.copy2(os.path.join(cfg.scriptdir, "rgetcmb"), ".") #may be okay as is


        #rstep3 is a move ( mv i* ../../lib_calib/<YYYYMM>/ )
        # and then a manual copy as the hetdex user to:
        # copy ../../lib_calib/<YYYYMM>/*     to   /scratch/projects/hetdex/lib_calib/<YYYYMM>


        #for rstep4
        shutil.copy2(os.path.join(cfg.scriptdir, "rgetcal1"), ".")

        #for rstep5 (note in the original method, rstep5 is automatically appended to rstep4)
        shutil.copy2(os.path.join(cfg.scriptdir, "rgetcal2"), ".")

        #rstep6 and rstep7 and rstep8 is a move and several copies
        # mv i* ../../lib_calib/<YYYYMM>/
        # as hetdex copy ../../lib_calib/<YYYYMM>/*     to   /scratch/projects/hetdex/lib_calib/<YYYYMM>
        # as hetdex copy ../../lib_calib/<YYYYMM>/*     to   /corral/utexas/Hobby-Eberly-Telesco/lib_calib/<YYYYMM>


        #shutil.copy2(os.path.join(cfg.scriptdir, "science_reductions", "rsetups"),".") #no, this function is its equivalent

        # shutil.copy2(os.path.join(cfg.scriptdir, "rfixspec"), ".")
        # shutil.copytree(os.path.join(cfg.scriptdir, "sciscripts"), ".", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir,"vdrp"), "vdrp", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir, "detect"), "detect", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir, "getcen"), "getcen", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir, "alldet"), "alldet", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir, "cs"), "cs", dirs_exist_ok=True)
        # shutil.copytree(os.path.join(cfg.scriptdir, "Diagnose"), "Diagnose", dirs_exist_ok=True)
        #
        # #update the "home" path tilde
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits")
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_fix")  # use '#' as sed separator rather than "/"
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_field")  # use '#' as sed separator rather than "/"
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_fix")  # use '#' as sed separator rather than "/"
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbacks")  # use '#' as sed separator rather than "/"
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_s")  # use '#' as sed separator rather than "/"
        # system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rimarb")  # use '#' as sed separator rather than "/"
        #
        #
        # #update the old red1 path; should now all point to the common directory above each sci<YYYYMMDDvSSS> directories
        # #yes, I want cwd_orig ... I want a single "reductions" directory off the top with all the sciXXX as siblings
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_field")  # use '#' as sed separator rather than "/"
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_fix")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rerun2")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rtaremc") #not necessary, just for completeness
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtar")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtarm.defunct") #not necessary, just for completeness
        # #next two are just safeties, the scripts THIS code is using should already have them updated
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# alldet/rfft")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# getcen/rgetifucen")
        #
        # #extra files needed
        # #vdrp : /work/00115/gebhardt/maverick/fplane ... need the fplane for the date
        # os.makedirs(os.path.join(cfg.cwd,"vdrp/fplane"), exist_ok=True)
        #
        #
        # if os.path.exists(os.path.join(karlfplane, f"fp{cfg.ifu_fqid[0:8]}")):
        #     shutil.copy2(os.path.join(karlfplane, f"fp{cfg.ifu_fqid[0:8]}"), os.path.join(cfg.cwd,"vdrp/fplane"))
        # else:
        #     print(f"{[cfg.ifu_fqid]} !Warning! fplane file (fp{cfg.ifu_fqid[0:8]}) not found. Using last known ({LastKnownFplane}) instead. ")
        #     shutil.copy2(os.path.join(karlfplane, f"{LastKnownFplane}"), os.path.join(cfg.cwd, f"vdrp/fplane/fp{cfg.ifu_fqid[0:8]}"))
        #
        # #fix paths in the . cfg files
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.original")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.gaia")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.sdss")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd}# vdrp/vdrp.config.panstarrs")
        #
        # system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config")
        # system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.original")
        # system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.gaia")
        # system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.sdss")
        # system_command(cfg, f"sed -i s#/work/03261/polonius/hetdex/science/sciscripts#{cfg.cwd}# vdrp/vdrp.config.panstarrs")
        #
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.original")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.gaia")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.sdss")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/vdrp.config.panstarrs")
        #
        # system_command(cfg, f"sed -i s#/scratch/00115/gebhardt#{cfg.cwd}# vdrp/shifts/runsh1")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/science_reductions#{cfg.cwd}# vdrp/shifts/runsh2")
        # system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# vdrp/shifts/run_shifts.sh")
        #
        # #vdrp: need the runshifts for the shot (from karlgettar)
        # #e.g. run_shifts.sh 20240730 009 16.317927 33.689304 1
        # # (formerly, this was the "clean_rta" script
        #
        # if os.path.exists(os.path.join(karlgettar,f"rta.{cfg.ifu_fqid[0:6]}")):
        #    rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(os.path.join(karlgettar,f"rta.{cfg.ifu_fqid[0:6]}"),
        #                                                         usecols=[1,2,3,4,5],unpack=True,dtype=str)
        # else:
        #     date_year = int(cfg.ifu_fqid[0:4])
        #     #if 202407 < int(cfg.ifu_fqid[0:6]) <= 202412:
        #     if date_year == 2024:
        #         rtafn = os.path.join(karlgettar,f"rta.202488")
        #     elif date_year >= 2025:
        #         rtafn = os.path.join(karlgettar, f"rta.{date_year}00")
        #         #rtafn = os.path.join(karlgettar, f"rta.202500")
        #     else:
        #         #should not happen
        #         rtafn = None #just to turn off warning
        #         Quit(cfg,-1,"Fatal. Cannot locate suitable rta file")
        #
        #     rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(rtafn,
        #                                                         usecols=[1, 2, 3, 4, 5], unpack=True, dtype=str)
        #
        # sel = (rta_date == cfg.ifu_fqid[0:8]) * (rta_shot == cfg.ifu_fqid[-3:])
        # with open(os.path.join(cfg.cwd,f"vdrp/shifts/rta.{cfg.ifu_fqid[0:6]}"),"w") as f:
        #     for d,s,ra,dec,v in zip(rta_date[sel], rta_shot[sel], rta_ra[sel], rta_dec[sel], rta_v[sel]):
        #         f.write(f"run_shifts.sh {d} {s} {ra} {dec} {v} \n")
        #

        #fix detect stuff

        return 0
    #end if not resume

    return 0


def node_setup(cfg): #,safelimit=0):
    """



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
            print(f"[{cfg.log_id}] checking if safe to start ...")

            while redlight:
                with lock:
                    #how many are active?
                    fns = glob.glob(os.path.join(Node_basedir, "*.sync"))
                    active = len(fns)
                    if active < MaxSafeActiveShots: #good to go
                        os.makedirs(Node_basedir, exist_ok=True)
                        with open(os.path.join(Node_basedir, f"{cfg.ifu_fqid}.sync"), "w") as f:
                            f.write(f"BEGIN {str(datetime.now())}\n")
                        redlight = False
                        print(f"[{cfg.log_id}] cleared to start.")
                    else:
                        print(f"[{cfg.log_id}] too many active shots ({active}). Must wait ...")

                if redlight:
                    time.sleep(SafeActiveShotsSleep)

        else:
            with lock:
                #create a sync directory and file
                #a bit later this will then be the count of simulataneous datevshots and used to tune the num of processes

                os.makedirs(Node_basedir, exist_ok=True)
                with open(os.path.join(Node_basedir,f"{cfg.ifu_fqid}.sync"), "w") as f:
                    f.write(f"BEGIN {str(datetime.now())}\n")

        # lock auto releases
    except:
        print(f"Exception! in node_setup()",traceback.format_exc())


def node_clean(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        lock = FileLock(Lock_tmp_mutex_fn)
        with lock:
            #clean up my sync
            try:
                #this might have already been removed (e.g. if post_clean() ran successfully)
                os.remove(os.path.join(Node_basedir,f"{cfg.ifu_fqid}.sync"))
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
            print(f"[{cfg.log_id}] checking active shot count for shared /tmp: {ct}")
        #lock auto releases
    except:
        print(f"Exception! in node_active_ct()", traceback.format_exc())

    return ct



def Quit(cfg,rc,msg=None,write_status=True):
    """

    :param cfg:
    :param rc:
    :param msg:
    :return:
    """

    if msg is not None:
        print(f"[{cfg.log_id}] ({rc})",msg)
    else:
        print(f"[{cfg.log_id}] ({rc})")


    if rc >=0:
        write_summary(cfg)

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
        system_command(cfg,f"rm {os.path.join(cfg.cwd,'status.run')}")
    except:
        pass

    try:
        if write_status and safe_cd(cfg.cwd):
            if rc < 0:
                with open("status.fail","w") as f:
                    f.write(f"[{cfg.log_id}] ({rc}) {msg}\n")
            else:
                with open("status.pass","w") as f:
                    f.write(f"[{cfg.log_id}] ({rc}) {msg}\n")
    except:
        pass

    node_clean(cfg)

    if rc >= 0:
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
            print(f"[{cfg.log_id}] subproc CMD: ({os.getcwd()}) > {cmd}")

        rc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

        if EchoCmds:
            print(f"[{cfg.log_id}] subproc rc = {rc}: CMD = ({os.getcwd()}) > {cmd}")
    except:
        print(f"[{cfg.log_id}] Exception! in blocking_command. cmd = {cmd}\n",traceback.format_exc())
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
        print(f"[{cfg.log_id}] CMD: ({os.getcwd()}) > {cmd}")

    if cfg.file_stdout:
        os.system(f"{cmd} &>> {cfg.file_stdout.name}")
    else:
        os.system(f"{cmd}")

def update_only(cfg):
    """

    :param cfg:
    :return:
    """

    print("*** Todo: fix paths for update_only()")

    lock = FileLock(Lock_mutex_fn)  # we are in the top directory (not sciXXXX)
    with lock:
        if cfg.update_local_repo:
            print("(Only) Updating local repo ...")
            shutil.copytree(os.path.join(ScriptRepo, ""),
                            os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
            cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)



def progress_update(cfg,progress_dict,key=None,status=True):
    """
    update the progress record
    :param progress_dict:
    :return:
    """

    try:
        if key is not None:
            progress_dict[key] = status

        fn = os.path.join(cfg.cwd,"progress.dat")
        with open(fn, 'w') as f:
            json.dump(progress_dict, f, indent=4)  # indent=4 for pretty-printing

    except:
        print(f"Exception in progress_update(). {traceback.format_exc()}")

def progress_init(cfg):
    """
    read/initialize progress dict
    :param cfg:
    :return:
    """
    dtprog = None
    try:
        fn = os.path.join(cfg.cwd,"progress.dat")
        if os.path.exists(fn):
            with open(fn, 'r') as f:
                dtprog = json.load(f)

            #assume this was successful, but also check --resume
            if cfg.resume is False:
                print("********** NOTICE **********")
                print("Call did NOT specfiy --resume, but this calibration appears to be at least partially complete.")
                print("Will attempt to resume (implied) at the last incomplete step ...")
                print("****************************")
                cfg.resume = True

    except:
        print(f"Exception in progress_init(). {traceback.format_exc()}")

    if dtprog is None:
        dtprog = {"rstep1": False,
                  "rstep2": False,
                  "rstep3": False, #manual move /scratch/projects
                  "rstep4": False,
                  "rstep5": False, #used to be part of rstep4
                  #"rstep6": False, #move to local lib_calib
                  #"rstep7": False, #manual move to /scratch/projects
                  #"rstep8": False, #manual move to /corral/
             }

        progress_update(cfg,dtprog) #write out the file, no updates yet

    return dtprog



def write_summary(cfg):
    """

    :param cfg:
    :return:
    """
    print("*******************")
    print("todo: write_summary")
    print("*******************")


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

        #always try to clean up /tmp
        node_clean(cfg)

        if cfg.clean <=0:
            print(f"[{cfg.log_id}] No -clean")
            cfg.clean_done = True
            return
        else:
            print(f"[{cfg.log_id}] --clean {cfg.clean}")


        if cfg.clean_only:  #this is at the top and cfg.cwd has NOT been set to cfg.ifu_fqid
            if not safe_cd(os.path.join(cfg.cwd,f"sci{cfg.ifu_fqid}")):
                print(f"[{cfg.log_id}] Could not initiate --clean")
                return
            else:
                #we did successfully jump to sciXXXX direcotry so update to that for the remainder
                cfg.cwd = os.getcwd()
        else:
            if not safe_cd(cfg.cwd):
                print(f"[{cfg.log_id}] Could not initiate --clean")
                return

        #os.chdir(os.path.join(cfg.cwd,f"sci{cfg.ifu_fqid}"))
        cfg.cwd = os.getcwd()  # now under the sci<shot> directory

        print(f"**************************")
        print(f"*** TODO  post_clean() ***")
        print(f"**************************")

        cfg.clean_done = True

    except:
        print(f"Exception in post_clean(). {traceback.format_exc()}")



def get_month_shot_info(cfg):
    """
    get the list of qualified shots for this YYYYMM

    these come from the gettar directory
    the rt*.YYYYMM files (or 202500, etc)
    and the runt<YYYYMM> files

    :param cfg:
    :return:
    """

    shotnum = []
    day = []
    try:

        fnbase, rt, run = get_calibration_tarlist(cfg)

        cfg.rtx_filelist  = rt
        cfg.runx_filelist = run
        cfg.yyyymm_fnbase = fnbase

        #need rt1.YYYYMM
        rt1 = os.path.join(karlgettar,f"rt1.{fnbase}")
        if rt1 in rt:
            with open(rt1,"r") as f:
                #example line:  run1t 20231018 003 exp01 202310
                for line in f:
                    if f"run1t {cfg.yyyymm}" in line:
                        shotnum.append(line.split()[2])
                        day.append(line.split()[1])

        else: #problem
            print(f"[{cfg.log_id}] rstep1 fail. Could not find {rt1}")
    except:
        print(f"Exception in get_month_shot_info(). {traceback.format_exc()}")

    return shotnum, day



def fetch_virus_tar(cfg):
    """

    copy from HETRaw_archive the <date>.tarfiles that correspond to what we need, untar and clean up

    :param cfg:
    :return:
    """

    rc = 0
    cwd = os.getcwd()
    try:

        one_at_a_time = True
        cmds = []
        log_stmts = []
        done_checks = []

        #example, extract specific file to specific location
        #tar -xf /corral/utexas/Hobby-Eberly-Telesco/het_raw/20230914.tar  20230914/virus/virus0000013.tar  -C ./20230914/virus/

        if safe_cd(f"{cfg.cwd_orig}/het_raw"):
            #copy each of them

            lockfn = f"mutex_{cfg.yyyymm}.lock"
            lock = FileLock(lockfn)
            with lock:
                #build up a list of commands to execute (optional copies, untars and a touch do indicate it is done)
                current_date = "19700101"
                cmd = ""
                for i in range(len(cfg.shotnum)):
                    #these are in date order
                    #if the files we expect are already there for this shot, skip and move on to the next
                    #expect to find <path>/<date>/virus/virus0000{shot}
                    expected_file = f"{cfg.cwd_orig}/het_raw/{cfg.day[i]}/virus/virus{str(cfg.shotnum[i]).zfill(7)}.tar"
                    if os.path.exists(expected_file):
                        print(f"[{cfg.log_id}] already copied file {i+1} of {len(cfg.shotnum)}: {expected_file}")
                        continue

                    #otherwise, extract the one file we want
                    fn = os.path.join(HETRaw_archive, f"{cfg.day[i]}.tar")
                    subtar_path = f"{cfg.day[i]}/virus/virus{str(cfg.shotnum[i]).zfill(7)}.tar"
                    output_path = f"./{cfg.day[i]}/virus/"

                    if one_at_a_time:
                        cmds.append(f"tar -xf {fn} {subtar_path} -C {output_path}")
                        log_stmts.append(f"[{cfg.log_id}] untarring {i + 1} of {len(cfg.shotnum)}: {fn}/{subtar_path} ...")
                    else: #all at once in batch
                        #if this is a new date OR the last iteration, write out the current command
                        if current_date != cfg.day[i]:
                            if len(cmd) != 0: #append the current command as a background task
                                cmd += " &"
                                cmds.append(cmd)
                                cmd = ""
                            #otherwise this is just the first comamnd
                            #start a new command
                            current_date = cfg.day[i]
                            cmd = f"tar -xf {fn} {subtar_path} -C {output_path} && " \
                                   f"touch {os.path.basename(fn)}_virus{str(cfg.shotnum[i]).zfill(7)}.done"
                            done_checks.append(f"{os.path.basename(fn)}_virus{str(cfg.shotnum[i]).zfill(7)}.done")
                        else:
                            cmd += f" ; tar -xf {fn} {subtar_path} -C {output_path} && " \
                                   f"touch {os.path.basename(fn)}_virus{str(cfg.shotnum[i]).zfill(7)}.done"
                            done_checks.append(f"{os.path.basename(fn)}_virus{str(cfg.shotnum[i]).zfill(7)}.done")

                #add the last command string chain, if needed
                if not one_at_a_time and len(cmd) > 0:
                    cmd += " &"
                    cmds.append(cmd)


                #build out the shell script and run
                if len(cmds) > 0:
                    if one_at_a_time:
                        print(f"[{cfg.log_id}] untar ({len(cmds)}), one at a time ... ",flush=True)
                        for log_stmt, cmd in zip(log_stmts,cmds):
                            print(log_stmt,flush=True)
                            system_command(cfg,cmd)

                    else: #all at once
                        with open(f"r{cfg.yyyymm}", "w") as f:
                            for cmd in cmds:
                                f.write(f"{cmd}\n")

                        done = False
                        ct = 10
                        while not done:
                            if os.path.exists(f"./r{cfg.yyyymm}"):
                                done = True
                            else:
                                ct -= 1
                                if ct > 0:
                                    time.sleep(1.0)
                                else:
                                    # this is a problem
                                    print(f"[{cfg.log_id}] Fatal. Cannot find tar fetch script r{cfg.yyyymm}",flush=True)
                                    done = True
                                    rc = -1

                        if rc >= 0:
                            system_command(cfg, f"chmod 755 r{cfg.yyyymm}")
                            system_command(cfg, f"./r{cfg.yyyymm}")
                        else:
                            return rc

                        done = False
                        while not done:
                            time.sleep(30.0)
                            for done_fn in done_checks:
                                if not os.path.exists(done_fn):
                                    break
                                done = True  # all accounted

                        # now clean up
                        for done_fn in done_checks:
                            try:
                                os.remove(done_fn)
                            except:
                                print(f"[{cfg.log_id}] Could not remove {done_fn}",flush=True)

                        try:
                            os.remove(f"r{cfg.yyyymm}")
                        except:
                            print(f"[{cfg.log_id}] Could not remove {cfg.yyyymm}",flush=True)

                else:
                    print(f"[{cfg.log_id}] No extra untarring to run.",flush=True)

            #end lock
        else:
            #this is a problem
            print(f"[{cfg.log_id}] temp het_raw location does not exist. Fatal.",flush=True)
            rc = -1

    except:
        print(f"Exception in fetch_virus_tar(). {traceback.format_exc()}")

    os.chdir(cwd)
    return rc

def rstep1(cfg):
    """

    :param cfg:
    :return:
    """

    rc = 0
    try:
        # example run1t 20231002 002 exp01 202310
        #         run1t <YYYYMMDD> <shotid> exp01 <YYYYMM>

        shotnum, day = get_month_shot_info(cfg)
        cfg.shotnum = shotnum
        cfg.day = day
        if len(shotnum) == 0: #this is fatal
            print(f"[{cfg.log_id}] rstep1  Could not get shotnumbers. Fatal")
            rc = -1

        #now set up the run1t calls
        #this would have been the 3 per line rt1.YYYYMM_1.run
        runfile = f"rt1.{cfg.yyyymm}_1.run"
        with open(runfile,"w") as f:
            for s,d in zip(shotnum,day):
                f.write(f"run1t {d} {s} exp01 {cfg.yyyymm} \n")


        #build the runt<YYYYMM> input
        with open(f"runt{cfg.yyyymm}","w") as f:
            ct = 0
            tot = len(day)
            for s, d in zip(shotnum, day):
                ct += 1
                f.write(f"echo {ct} / {tot} ; rback {d} {s} exp01 {cfg.ifuslot} {cfg.specid} {cfg.yyyymm} 0 \n")

        #get the shot tar files and extract them
        rc = fetch_virus_tar(cfg)

        #now, we run it
        if rc >= 0:
            system_command(cfg, f"chmod 755 {runfile}")
            system_command(cfg,f"./{runfile}")


    except:
        print(f"Exception in rstep1(). {traceback.format_exc()}")

    return rc



def rstep2(cfg):
    """

    :param cfg:
    :return:
    """

    rc = 0
    try:
        #just one call (only 1 IFU)
        system_command(cfg,f"rgetcal0 {cfg.ifuslot} {cfg.yyyymm} {cfg.specid}")
    except:
        print(f"Exception in rstep2(). {traceback.format_exc()}")

    return rc


def rstep3(cfg):
    """

    two parts
    1)   mv i* ../../lib_calib/<YYYYMM>/.    (as local user)
    2) also copy  to /scratch/projects/hetdex/lib_calib/<YYYYMM>  as hetdex user


    :param cfg:
    :return:
    """

    rc = 0
    try:

        #print(f"rstep3, mv i* to local lib_callib")

        #we only want certain files ... and lets copy instead of move
        #
        #system_command(cfg,f"mv i* ../../lib_calib/{cfg.yyyymm}")

        system_command(cfg, f"cp i*cbwt* ../../lib_calib/{cfg.yyyymm}")
        system_command(cfg, f"cp i*cbxt* ../../lib_calib/{cfg.yyyymm}")
        system_command(cfg, f"cp i*cbmf* ../../lib_calib/{cfg.yyyymm}")
        system_command(cfg, f"cp i*cbmp* ../../lib_calib/{cfg.yyyymm}")
        system_command(cfg, f"cp i*wave* ../../lib_calib/{cfg.yyyymm}")

    except:
        print(f"Exception in rstep3(). {traceback.format_exc()}")

    return rc


def rstep4(cfg):
    """

    much like rstep1, same format of calls but for rgetcal1 instead of run1t

    :param cfg:
    :return:
    """

    rc = 0
    try:
        # example run1t 20231002 002 exp01 202310
        #         run1t <YYYYMMDD> <shotid> exp01 <YYYYMM>

        shotnum, day = get_month_shot_info(cfg)
        cfg.shotnum = shotnum
        cfg.day = day
        if len(shotnum) == 0: #this is fatal
            print(f"[{cfg.log_id}] rstep4  Could not get shotnumbers. Fatal")
            rc = -1

        #now set up the run1t calls
        #this would have been the 2 per line rt2.YYYYMM_1.run
        runfile = f"rt2.{cfg.yyyymm}_1.run"
        with open(runfile,"w") as f:
            for s,d in zip(shotnum,day):
                f.write(f"rgetcal1 {d} {s} exp01 {cfg.yyyymm} 1\n")
                #last charater is 0 (new calibration) or 1 (read in a calibration)
                #in this case, we want the 1 (the calibration bit was done in prior steps)

        Quit(cfg,99,"*** TESTING *** \n early exit in rstep4, after create rt2.xxx_1.run")

        #now, we run it
        if rc >= 0:
            system_command(cfg, f"chmod 755 {runfile}")
            system_command(cfg,f"./{runfile}")


    except:
        print(f"Exception in rstep1(). {traceback.format_exc()}")

    return rc

def rstep5(cfg):
    """

    in the old way, this was appended to rstep4, here it is split back out

    :param cfg:
    :return:
    """

    rc = 0
    try:
        #just one call (only 1 IFU)
        system_command(cfg,f"rgetcal2 {cfg.ifuslot} {cfg.yyyymm} {cfg.specid}")
    except:
        print(f"Exception in rstep5(). {traceback.format_exc()}")

    return rc

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

if cfg.clean_only:
    print(f"[{cfg.log_id}] Performing only the CLEAN, level : {cfg.clean} ...")
    post_clean(cfg)
    Quit(cfg,0,"Clean complete. Exiting",write_status=False)

rc = precheck(cfg)
if rc < 0:
    Quit(cfg,rc,"FATAL! Precheck failed. Reduction cannot run.",write_status=False)


rc = initial_setup(cfg)

if rc < 0:
    Quit(cfg,rc,"Could not complete initial setup.",write_status=False)

print("todo: node_setup")
#rc = node_setup(cfg)

if cfg.special == 1:
    print("*** TODO Special Handling? ***")
    Quit(cfg, 0, f"Done with special handling. {cfg.ifu_fqid}",write_status=False)


#########
# after the initial setup, move stdout and stderr to a log file
#########

cfg.file_stdout = open(f"{cfg.ifu_fqid}.log","a")
print(f"[{cfg.log_id}] Logging redirected to: {cfg.cwd}/{cfg.file_stdout.name}")


# get the progress state. Useful if resuming (implied)
dtprog = progress_init(cfg)


#begin
with open("status.run", "w") as f:
    f.write(f"[{cfg.log_id}] running .... \n")


#todo: for rstep1 on vm-small can run 4 maybe 5 at a time. Around 6-8 GB RAM needed per
# so could use mutex counting (like on reduce_shot.py), but need to be careful on the rt1.xxx_1.run calls since
# those are in a shell script ... would have to change s|t each line in rt1.xxx_1.run is called from Python and
# not from the rt1.xx_1.run file

if not dtprog["rstep1"]:
    rc = rstep1(cfg)
    if rc < 0:
        Quit(cfg,rc,"Could not complete rstep1.",write_status=False)
    else:
        progress_update(cfg, dtprog, "rstep1")
else:
    print(f"[{cfg.log_id}] Skipping rstep1")

if not dtprog["rstep2"]:
    rc = rstep2(cfg)
    if rc < 0:
        Quit(cfg, rc, "Could not complete rstep2.", write_status=False)
    else:
        progress_update(cfg, dtprog, "rstep2")
else:
    print(f"[{cfg.log_id}] Skipping rstep2")

if not dtprog["rstep3"]:
    rc = rstep3(cfg)
    if rc < 0:
        Quit(cfg, rc, "Could not complete rstep3.", write_status=False)
    else:
        progress_update(cfg, dtprog, "rstep3")

        print(f"*** User ***")
        print(f"You now need to change to the 'hetdex' user and copy the contents of:  ")
        print(f"    {cfg.cwd_orig}/lib_calib/{cfg.yyyymm} ")
        print(f"    to")
        print(f"    /scratch/projects/hetdex/lib_calib/{cfg.yyyymm}")
        print(f"    and")
        print(f"    /corral/utexas/Hobby-Eberly-Telesco/lib_calib/{cfg.yyyymm}")
        print(f"\n")
        print(f"This update to the calibration is done. The rest of the steps only apply for a FULL IFU array calibration.")
        #print(f"You can then re-run calibrate_ifu.py BUT you must add --resume ")
        #
        #Quit(cfg, 0, f"[{cfg.log_id}] Partially Complete. Waiting user action.")
else:
    print(f"[{cfg.log_id}] Skipping rstep3")



#
# !!! we may actually be done ... the next steps combine and check across ALL the IFUs
#     (e.g. Amp to Amp (ata) stuff), which we CANNOT do since this is only for one IFU
#

#
# if not dtprog["rstep4"]:
#     rc = rstep4(cfg)
#
#     if rc < 0:
#         Quit(cfg, rc, "Could not complete rstep4.", write_status=False)
#     else:
#         progress_update(cfg, dtprog, "rstep4")
# else:
#     print(f"[{cfg.log_id}] Skipping rstep4")
#
#
# if not dtprog["rstep5"]:
#     rc = rstep5(cfg)
#
#     if rc < 0:
#         Quit(cfg, rc, "Could not complete rstep5.", write_status=False)
#     else:
#         progress_update(cfg, dtprog, "rstep5")
# else:
#     print(f"[{cfg.log_id}] Skipping rstep5")

##########
# DONE
##########

post_clean(cfg)

Quit(cfg, 0, f"[{cfg.log_id}] Complete")