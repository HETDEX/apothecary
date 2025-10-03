"""
Perform science reduction (astrometry, flux calibration, line and continuum detection, elixer )
  on a single shot (field). Can be HETDEX dither-style or not.

  Copies code source and works in the current directory

  input: [required] shotid (or datevshot)
         [optional] -clean  (clean up workfiles, leaving only the output and logs)
         [optional] -overwrite (delete and overwrite the datevshot directory)
         [optional] -exp <##> (specify a single exposure to reduce, if there is more than one in the shot; 0 = all)



This is just a large Python script. There is no defined main(), but the principle logic begins in a
  "main" commented section.

Error control (at least for now) is deliberately limited as I want no hidden errors. Mostly if anything is wrong
   I want this to break.

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
from h5tools import amp_stats as AmpStats
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

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
plt.style.use('default')

########################################################################
# CONFIGURATION
########################################################################
EchoCmds = True  #if True echo system commands to the log
FilterDetsOnBadAmps = False # if True, do NOT pass detections that are on reported bad amps to elixer for processing
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

HETRaw_archive = "/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"
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




#execute steps
s01_run1s = True

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

    clean: int = 0 #post run clean level; 0 = do not clean
    clean_only: bool = False #if True, do nothing except for -clean
    #simul : int = 0 #number of simultaneous shots being run (e.g. tasks per node), 0 is unset
    update_local_repo: bool = False
    overwrite: bool = False
    resume: bool = False
    shotid: int = 0
    datevshot: str = ""
    exp: int = 0  #specific exposure number to reduce
    numexp: int = 0 #number of exposures in the shot
    cwd_orig: str = os.getcwd()
    cwd: str = os.getcwd()
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




########################################################################
# Basic user input
########################################################################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

cfg = Config()

if "-help" in args:
    #do NOTHING else except print the help
    help = """
    
    usage: python reduce_shot.py <shot> [switches]
           where <shot> is YYYYMMDDSSS or YYYYMMDDvSSS
    
    switches (all optional):
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
                      
    --exp <integer> : will operate on only the specified exposure (e.g. in a multi-exposure observation, can select
                      exactly one to reduce). If not present or set to (0), will use all exposures for the observation. 
        
    --help : display this help text and exit
             
    --overwrite : removes the shot working directory completely and (re)starts fresh. 
              !!! Notice: --resume has priority over --overwrite
    
    --resume : (re)starts roughly at the last completed step (see sciXXXX/progress.dat)
           !!! Notice: This does NOT re-run steps that completed with failures, it only re-runs incomplete steps.
           !!! Notice: --resume has priority over --overwrite
               
    --update : removes and re-fetches the local_script_depo prior to running
               on a --resume, also updates the scripts already in the shot working directory
    
    """

    print(help)

    exit(0)

queue_elixer = False
if "-queue_elixer" in args:
    print("Hidden switch : queueing elixer slurm jobs that match datevshot ...")
    args.remove("-queue_elixer")
    queue_elixer = True

if "-clean" in args:
    i = args.index("-clean")
    try:
        cfg.clean = int(args[i+1])
        if cfg.clean < 0:   #negative values are the same, but force just the clean operation (e.g. to be used on an old run)
            cfg.clean *= -1
            cfg.clean_only = True
    except:
        print(f"Invalid -clean specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-clean")
else:
    cfg.clean = DefaultClean #usually this is level 1


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


#whatever is left should be the shot
if len(args) != 1:
    print(f"Fatal: Problem with remaining args: {args}")
    print(f"exititing....")
    exit(-1)
else:
    if not queue_elixer:
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
    #print(fns)
    for fn in fns:
        try:
            for dettype in ["line","cont"]:
                slurm_path = os.path.join(fn, f"elixer/{dettype}")
                if safe_cd(slurm_path):
                    if os.path.exists("elixer.slurm"):
                        cmd = f"cd {slurm_path} ; sbatch elixer.slurm ; cd {cwd}"
                        #print(cmd)
                        system_command(cfg,cmd)

                    os.chdir(cwd)
        except:
            print(traceback.format_exc())

    os.chdir(cwd)

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
        #always try to clean up /tmp
        node_clean(cfg)

        if cfg.clean <=0:
            print(f"[{cfg.datevshot}] No -clean")
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
                    "rback_fix","rbacks","rbfits","rbfits0","rbfits_fix","rbfits_s","rerun","rerun2","rfield",
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
            # getcen
            ############################
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf ./getcen")


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
                system_command(cfg, f"rm 20250828011_stats.pickle")
                system_command(cfg, f"rm 20250828011_ampstats.fits")

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
                ############################
                # individual exposures
                ############################
                system_command(cfg, f"rm -rf d{cfg.datevshot.replace('v', 's')}exp*")

                system_command(cfg, f"rm -rf alldet")
                system_command(cfg, f"rm -rf cs")



        if cfg.clean >= 5:  # ONLY keep the shot
            if safe_cd(cfg.cwd):
                system_command(cfg, f"rm -rf elixer") #NOTE: here we can go ahead and remove since it is the top dir
                                                      #and does not matter if elixer has executed
                system_command(cfg, f"rm diagnose_classifications.tab")
                system_command(cfg, f"rm {cfg.datevshot}_*")


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
                print("Call did NOT specfiy --resume, but this reduction appears to be at least partially complete.")
                print("Will attempt to resume (implied) at the last incomplete step ...")
                print("****************************")
                cfg.resume = True

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

        progress_update(cfg,dtprog) #write out the file, no updates yet

    return dtprog

def Quit(cfg,rc,msg=None,write_status=True):
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

    try:
        if write_status and safe_cd(cfg.cwd):
            if rc < 0:
                with open("status.fail","w") as f:
                    f.write(f"[{cfg.datevshot}] ({rc}) {msg}")
            else:
                with open("status.pass","w") as f:
                    f.write(f"[{cfg.datevshot}] ({rc}) {msg}")
    except:
        pass

    node_clean(cfg)

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


def get_guider_fwhm(cfg):
    """
    try to get the seeing fwhm (IQ -- image quality) from the guider (gc1 or gc2) for the shot
    :param cfg:
    :return:
    """

    try:
        saved_fn = "guider.fwhm"
        if os.path.exists(saved_fn):
            fwhm = np.loadtxt(saved_fn,dtype=float) #just one value
            if fwhm is None or np.isnan(fwhm):
                 #we will continue below and rebuild
                print(f"Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
            else:
                fail = False
                try:
                    if len(fwhm) == 0 or fwhm[0] <= 0:
                        print(f"Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
                        fail = True
                except:
                    pass

                try:
                    fwhm = float(fwhm)
                except:
                    print(f"Found {saved_fn}, but guider seeing FWHM is inavlid ({fwhm}). Will recompute ...")
                    fail = True

                if not fail:
                    print(f"Found {saved_fn}. Using guider seeing FWHM = {fwhm}")
                    return fwhm

        exposure_times = [] #exposure times (seconds) from HDU
        exposure_fn_times =[] #date T time string from the fileanames

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
        base_tarfn = os.path.join(path, f"virus/{virus_shot}.tar")
        if os.path.exists(base_tarfn):
            with tar.open(base_tarfn, "r") as tarfh:
                fns = np.array(tarfh.getnames())
                #should look like a list of:  virus0000007/exp01/virus/20241017T024257.2_106RU_sci.fits
                #should all be the same base within each exposure
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
                        exposure_times.append(fh[0].header['EXPTIME'])
                        #or would DARKTIME be better ??
                        fh.close()


        if len(exposure_fn_times) == 0:
            #could not find any
            return None

        #try gc1
        base_tarfn = os.path.join(path,"gc1/gc1.tar")
        if os.path.exists(base_tarfn):
            with tar.open(base_tarfn, "r") as tarfh:
                gc1_names = tarfh.getnames()
                #should look like a list of:  20241017T024256.0_gc1_sci.fits
                #only want the dateTtime prefix
                #gc1_times = [x.split("/")[1].split("_")[0] for x in gc1_names]

        #try gc2
        base_tarfn = os.path.join(path,"gc2/gc2.tar")
        if os.path.exists(base_tarfn):
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

            for fntime, exptime in zip(exposure_fn_times, exposure_times):
                fntime = datetime(int(fntime[0:4]), int(fntime[4:6]), int(fntime[6:8]),
                                  int(fntime[9:11]), int(fntime[11:13]), int(fntime[13:15]))
                delta_time = fntime - timedelta(seconds=exptime)
                time_str = f"{str(delta_time.year)}{str(delta_time.month).zfill(2)}{str(delta_time.day).zfill(2)}T" \
                           f"{str(delta_time.hour).zfill(2)}{str(delta_time.minute).zfill(2)}{str(delta_time.second).zfill(2)}.0"

                gc_start_times.append(time_str)
                #Guider typically records every few seconds, so no need to go after
                delta_time = fntime #+ timedelta(seconds=120.0)  # go up to 2 minutes after
                time_str = f"{str(delta_time.year)}{str(delta_time.month).zfill(2)}{str(delta_time.day).zfill(2)}T" \
                           f"{str(delta_time.hour).zfill(2)}{str(delta_time.minute).zfill(2)}{str(delta_time.second).zfill(2)}.0"
                gc_stop_times.append(time_str)

            if gc1_names is not None:
                gc1_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                for start_time, stop_time in zip(gc_start_times,gc_stop_times):
                    for name in gc1_names:
                        #just want the time part
                        if start_time <= name[2:19] <= stop_time: #yes, technically string compares, but works for the time format here
                            gc1_near.append(name)

            if gc2_names is not None:
                gc2_near = [] #list of lists ... ie. if 3 exposures, will be a 3 long list each with 2 elements
                for start_time, stop_time in zip(gc_start_times,gc_stop_times):
                    for name in gc2_names:
                        #just want the time part
                        if start_time <= name[2:19] <= stop_time: #yes, technically string compares, but works for the time format here
                            gc2_near.append(name)
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
        if GuiderFWHM_ALL:
            #note: which guider is active could change? so check both?
            base_tarfn = os.path.join(path, "gc1/gc1.tar")
            if os.path.exists(base_tarfn):
                for name in gc1_near:
                    try:
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, name)
                        fh = fits.open(t1)
                        if fh[0].header['GUIDLOOP'] == 'ACTIVE':
                            iq.append(float(fh[0].header['IQ'])) #this is clipped later to exlcude 0s and other bad values
                        fh.close()
                        t1.close()
                    except:
                        pass #don't let one bad read bomb out

            base_tarfn = os.path.join(path, "gc2/gc2.tar")
            if os.path.exists(base_tarfn):
                for name in gc2_near:
                    try:
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, name)
                        fh = fits.open(t1)
                        if fh[0].header['GUIDLOOP'] == 'ACTIVE':
                            iq.append(float(fh[0].header['IQ'])) #this is clipped later to exlcude 0s and other bad values
                        fh.close()
                        t1.close()
                    except:
                        pass #don't let one bad read bomb out

        else: #just the two nearest
            #now which to use? gc1 or gc2
            base_tarfn = os.path.join(path, "gc1/gc1.tar")
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


            base_tarfn = os.path.join(path, "gc2/gc2.tar")
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
                base_tarfn = os.path.join(path, "gc1/gc1.tar")
                gc_near = gc1_near

            elif gc2_active:
                base_tarfn = os.path.join(path, "gc2/gc2.tar")
                gc_near = gc2_near
            else:
                print(f"No active guider. Cannot get seeing fwhm.")
                return None

            for near_list in gc_near:
                for near_file in near_list:
                    try:
                        t1, p1 = Utils.open_file_from_tar(base_tarfn, near_file)
                        fh = fits.open(t1)
                        iq.append(float(fh[0].header['IQ']))
                        fh.close()
                        t1.close()
                    except:
                        print(f"Invalid IQ card")

        if len(iq) == 1:
            np.savetxt(saved_fn,[iq[0]],fmt="%0.4f")
            return iq[0]
        elif len(iq) > 1:
            md_iq = np.nanmedian(np.clip(iq,0.1,9.0))
            np.savetxt(saved_fn, [md_iq], fmt="%0.4f")
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
        month = cfg.datevshot[0:6]
        path_check = os.path.join(hetdex_projects_path,f"lib_calib/{month}")
        if not os.path.exists(path_check):
            print(f"Precheck fail. Dir does not exist: {path_check}")
            return -1

        #common missing installs (that don't show up until later)
        pkgs=['sklearn','numba']
        rc = 0
        for pkg in pkgs:
            if importlib.util.find_spec(pkg) is None:
                print(f"Fatal. You need to install '{pkg}'")
                rc = -1

        if rc != 0:
            return -1

    except:
        print(f"Exception in precheck: {traceback.format_exc()}")

    return 0

def initial_setup(cfg):
    """
    copy from script repo(s)

    change the cwd to the new workdir for this shot

    this is largely equivalent to what rsetups would do, but here for a single shot and not a whole month

    :return:
    """

    workdir = os.path.join(WorkDirRoot,f"sci{cfg.datevshot}")

    resume = False #notice: cfg.resume MAY be true and this can still be false if the directory does not already exist
    if os.path.exists(workdir):
        if cfg.resume:
            print(f"[{cfg.datevshot}] Resuming. Leave directory intact: {workdir}")
            resume = True
        elif cfg.overwrite:
            print(f"[{cfg.datevshot}] Overwriting directory {workdir} ... ")
            shutil.rmtree(workdir)
        else:
            print(f"[{cfg.datevshot}] Shot directory already exists here! {workdir}")
            print(f"[{cfg.datevshot}] Please include --resume or --overwrite to make intention clear.")
            print(f"[{cfg.datevshot}] Or is this a repeated SLURM task?")
            return -1

    if not resume:
        os.makedirs(workdir)

        if LocalScriptRepo is not None:
            #need to see if can get the file lock in case another instance is copying
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:
                if cfg.update_local_repo:
                    print("Updating local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(), LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

                if os.path.exists(LocalScriptRepo): #we want to use it
                    print("Using local repo ...")
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
                else:
                    #copy first to local script repo
                    print("Copying to local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                    os.path.join(os.getcwd(),LocalScriptRepo), dirs_exist_ok=True)
                    cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)

            #lock auto releases
        else:
            print("Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")
    else:
        if LocalScriptRepo is not None:
            fatal_rtn = False
            lock = FileLock(Lock_mutex_fn) #we are in the top directory (not sciXXXX)
            with lock:

                if cfg.update_local_repo:
                    print("Updating local repo ...")
                    shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
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
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")


    os.chdir(workdir)
    cfg.cwd = os.getcwd() #now under the sci<shot> directory


    #check for virus.tar files
    virustar = f"{cfg.datevshot[0:8]}/virus/virus0000{cfg.datevshot[-3:]}.tar"
    virus_paths = [HET_by_date,os.path.join(cfg.cwd_orig,"het_raw")]
    if not os.path.exists(os.path.join(virus_paths[0],virustar)):
        if not os.path.exists(os.path.join(virus_paths[1],virustar)):
            fail = True
            if os.path.exists(os.path.join(HETRaw_archive,f"{cfg.datevshot[-3:]}.tar")):
                print(f"[{cfg.datevshot}]. Could not locate {virustar} under {virus_paths}")
                print(f"[{cfg.datevshot}]. Will attempt to copy. If successful, this will take many minutes ...")

                #tmp lock file (so only one attempt to copy and extract; since there are often many observations
                #   in this file, only one of the tasks that are on that date should copy)
                tmp_lock_file =f"{cfg.datevshot[-3:]}.lock"
                lock = FileLock(tmp_lock_file)  # we are in the top directory (not sciXXXX)
                with lock:
                    #try to copy
                    if not os.path.exists(f"{cfg.cwd_orig}/het_raw"):
                        os.makedirs(f"{cfg.cwd_orig}/het_raw", exist_ok=True)

                    if safe_cd("../het_raw"):
                        #check that another process has not already copied the file (see above pathing)
                        if not os.path.exists(os.path.join(virus_paths[1], virustar)):

                            shutil.copy2(os.path.join(os.path.join(HETRaw_archive,f"{cfg.datevshot[-3:]}.tar"),"."))

                            try:
                                # need only some of the paths ... don't know which /virus we may also need so keep all
                                cmd = f"tar -xvf {cfg.datevshot[-3:]}.tar {cfg.datevshot[-3:]}/virus"
                                system_command(cfg, cmd)

                                cmd = f"tar -xvf {cfg.datevshot[-3:]}.tar gc1"
                                system_command(cfg, cmd)

                                cmd = f"tar -xvf {cfg.datevshot[-3:]}.tar gc2"
                                system_command(cfg, cmd)

                                #last check (note: gc1 and gc2 are desirable to have, but not required)
                                if os.path.exists(os.path.join(virus_paths[1], virustar)):
                                    fail = False #we should be okay now

                            except:
                                print(f"[{cfg.datevshot}] Failed to copy/extract date tar file.",traceback.format_exc())

                            #now we should delete the big <date>.tar file
                            cmd = f"rm  {cfg.datevshot[-3:]}.tar"
                            system_command(cfg, cmd)

                        else:
                            fail = False #another process DID copy and the file we want is there now

            if fail:
                print(f"[{cfg.datevshot}] FATAL. Could not locate {virustar} under {virus_paths}")
                print(f"[{cfg.datevshot}] You may need to first copy and extract "
                      f"{HETRaw_archive}/<date>.tar to your local het_raw directory")
                return -1
            else: #we did eventually get what we need, so go back to the working dir and continue
                os.chdir(workdir)

    if not resume or cfg.update_local_repo:
        print(f"Copying source code to working directory {cfg.cwd}...")
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
            print(f"{[cfg.datevshot]} !Warning! fplane file (fp{cfg.datevshot[0:8]}) not found. Using last known ({LastKnownFplane}) instead. ")
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
           rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(os.path.join(karlgettar,f"rta.{cfg.datevshot[0:6]}"),
                                                                usecols=[1,2,3,4,5],unpack=True,dtype=str)
        else:
            if 202407 < int(cfg.datevshot[0:6]) <= 202412:
                rtafn = os.path.join(karlgettar,f"rta.202488")
            elif int(cfg.datevshot[0:6]) >= 202500:
                rtafn = os.path.join(karlgettar, f"rta.202500")
            else:
                #should not happen
                rtafn = None #just to turn off warning
                Quit(cfg,-1,"Fatal. Cannot locate suitable rta file")

            rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(rtafn,
                                                                usecols=[1, 2, 3, 4, 5], unpack=True, dtype=str)

        sel = (rta_date == cfg.datevshot[0:8]) * (rta_shot == cfg.datevshot[-3:])
        with open(os.path.join(cfg.cwd,f"vdrp/shifts/rta.{cfg.datevshot[0:6]}"),"w") as f:
            for d,s,ra,dec,v in zip(rta_date[sel], rta_shot[sel], rta_ra[sel], rta_dec[sel], rta_v[sel]):
                f.write(f"run_shifts.sh {d} {s} {ra} {dec} {v} \n")


        #fix detect stuff
        system_command(cfg, f"sed -i s/rsetstar/#/ detect/rallcal") #comment out rsetstar call as we want to run it separately
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rfitfw0")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsetstar")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsp3f")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/single_shot/science_reductions#{cfg.cwd}# detect/rsp3fc")


        return 0
    #end if not resume

    return 0


def node_setup(cfg):
    """

    extra stuff shared for datevshots on same node

    :param cfg:
    :return:
    """

    try:
        lock = FileLock(Lock_tmp_mutex_fn)
        with lock:
            #create a sync directory and file to hold list of datevshot being used
            #a bit later this will then be the count of simulataneous datevshots and used to tune the num of processes

            os.makedirs(Node_basedir, exist_ok=True)
            with open(os.path.join(Node_basedir,f"{cfg.datevshot}.sync"), "w") as f:
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
                os.remove(os.path.join(Node_basedir,f"{cfg.datevshot}.sync"))
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
        lock = FileLock(Lock_tmp_mutex_fn)
        with lock:
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
        idx = -1
        date, shot, exp = np.loadtxt(fn.replace('?','s'),usecols=[1,2,3],unpack=True,dtype=str)
        dse = sorted(np.unique([d + s + e for d, s, e in zip(date, shot, exp)]))
        ds = np.array([x[:-5] for x in dse])
        ct = np.count_nonzero(ds == str(shotid))
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

    #otherwise we did find it ... how may exposures?
    return ct, fn


def run_run1s(cfg):
    """

    science_reduction/sciscripts/run1  batch file


    :return:
    """


    print("run1s ...")
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

        #cmd = "sed -i s#\${scriptdir}"+f"#{cfg.scriptdir}/sciscripts/# run1s"
        #scripts have already been copied to shot workding dir
        cmd = "sed -i s#\${scriptdir}" + f"#{cfg.cwd}/# run1s"
        system_command(cfg,cmd)  # use '#' as sed separator rather than "/"

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

    print("todo: check run1s ... see check_red1.ipynb")
    return 0


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

        fns = sorted(fns)
        for fn in fns:
            norms = np.loadtxt(os.path.join(fn, "norm.dat"))  # one line, 3 values
            if np.count_nonzero(abs(1 - norms) > 0.5) > 0 or np.any(np.isnan(norms)):
                print("Possible bad dither norm:", os.path.basename(fn), norms)
                rc = -1

        print("check norms ... fail")
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


    #change to the vdrp/shifts directory
    # since we do not know if we will need all of GAIA, SDSS and PanSTARRS at this point,
    #    run all three of them


    os.chdir(os.path.join(cfg.cwd,"vdrp/shifts"))


    #GAIA first
    #command is based on rta.YYYYMM
    #run_shifts 20240730 009 16.317927 33.689304 1  20240730 GAIA

    print(f"[{cfg.datevshot}] VDRP: GAIA")
    try:
        with open(f"rta.{cfg.datevshot[0:6]}","r") as rta:
            for line in rta: #really should only be one line
                if len(line) > 10:
                    #we DON'T want the carriage return
                    line = line.rstrip('\n')
                    system_command(cfg,f"{line} {cfg.datevshot[0:8]} GAIA")


        vdrp_check_norms(cfg)
        vdrp_check_shout_ifu(cfg)
        vdrp_cp2dithall(cfg,"gaia")
        #move the output under gaia
        os.makedirs("gaia",exist_ok=True)
        system_command(cfg,f"mv {cfg.datevshot[0:6]}??v??? gaia")
    except Exception as e:
        print(f"[{cfg.datevshot}] VDRP: GAIA fail.",e,  "\n", traceback.format_exc())


    #!!! notice: we are doing limited checking here ... that will be done later
    # previously would run: check_norms YYYYMM   (included)
    #                       check_shot.ifu YYYYMM (included)
    #  examine the .pngs manually  NOT DONE
    #                  run: make_good_shots YYYYMM  NOT DONE
    #                       cp2dithall YYYYMM  (included)


    #SDSS
    print(f"[{cfg.datevshot}] VDRP: SDSS")
    try:
        with open(f"rta.{cfg.datevshot[0:6]}","r") as rta:
            for line in rta: #really should only be one line
                if len(line) > 10:
                    #we DON'T want the carriage return
                    line = line.rstrip('\n')
                    system_command(cfg,f"{line} {cfg.datevshot[0:8]} SDSS")

        vdrp_check_norms(cfg)
        vdrp_check_shout_ifu(cfg)
        vdrp_cp2dithall(cfg,"sdss")
        #move the output under sdss
        os.makedirs("sdss",exist_ok=True)
        system_command(cfg,f"mv {cfg.datevshot[0:6]}??v??? sdss")
    except Exception as e:
        print(f"[{cfg.datevshot}] VDRP: SDSS fail.", e, "\n", traceback.format_exc())


    if do_panstarrs:
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
        except Exception as e:
            print(f"[{cfg.datevshot}] VDRP: PANSTARRS fail.", e, "\n", traceback.format_exc())

    #todo: which is the main dithall, etc???? (normally it is SDSS for calibration and gaia for astrometry)
    #print("!!! todo: copy GAIA dithaall to /scatch/projects and /corral-repl ???")

    # set the star catalog to use for calibration (make a softlink to the starcat specific output for this shot)
    # just in case something went badly wrong, make sure we are in the right directory
    os.chdir(os.path.join(cfg.cwd, "vdrp/shifts"))
    system_command(cfg, f"ln -s {cfg.starcat_cal}/{cfg.datevshot} {cfg.datevshot}")

    os.chdir(cfg.cwd)


def check_vdrp(cfg):
    """

    :param cfg:
    :return:
    """
    print("todo: check vdrp ... ")
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

    print(f"[{cfg.datevshot}] flux calibration ... ")

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

    print(f"[{cfg.datevshot}] flux calibration (rsetstar) ... ")
    #call rsetstar independently
    system_command(cfg, f"rsetstar {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {star_catalog}")

    print(f"[{cfg.datevshot}] flux calibration (rallcal) ... ")
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
            if badct> 0:
                rc = -1
                print(f"[{cfg.datevshot}] failed to get IFU centers. Invalid RA, Decs.")
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
                        print(f"[{cfg.datevshot}] rff [Pass] {fn}")
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
            base_str = f"{ct}) checking {multi[6:]}  {ra:0.7f} {dec:0.7f} ... "
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
            print(f"{cfg.datevshot} (run_rcal) All Pass")
        elif len(failed_rcal_list) == len(ras):  # all failed
            rc = -1
            print(f"{cfg.datevshot} (run_rcal) ALL FAIL")
        else:
            rc = 1
            print(f"{cfg.datevshot} (run_rcal) Mixed results of {len(ras)}: {len(passed_rcal_list)} Pass, {len(failed_rcal_list)} FAIL")


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
        cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 104\n"
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
        num_procs = max(1,min(num_procs, MaxPerShotProcs_mp_rf1))

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

        #now - check them all (this can be serial, it is fast)
        ct = 0
        for multi, ra, dec in zip(multis, ras, decs):
            ct +=1
            base_str = f"{ct}) checking [{cfg.datevshot}] {multi[6:]}  {ra:0.7f} {dec:0.7f} ... "
            #print(f"{ct}) checking [{cfg.datevshot}] {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ", end="")
            output_found = np.array([0,0,0])
            for i,ext in enumerate(output_extensions):
                outfn = f"{cfg.datevshot}_{multi[6:]}{ext}"

                if os.path.exists(os.path.join("detect_out/", outfn)):
                    output_found[i] = 1

            if np.count_nonzero(output_found) != 3:
                #something failed, we will want to re-run these once
                cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 104\n"
                failed_list.append(cmd)
                print(f"{base_str} FAIL. May re-run at the end.")
            else:
                print(f"{base_str} pass")

        #let any that failed re-run as serial to keep it simple
        if len(failed_list) > 0:
            if len(failed_list) < len(ras):
                #print(f"{base_str} {len(failed_list)} failed. Can be transient issues, so will re-run ...")
                print(f"{cfg.datevshot} (rdet_rf1) : serially re-run {len(failed_list)} failed IFUs ...")
                for ct,cmd in enumerate(failed_list):
                    print(f"{cfg.datevshot} (rdet_rf1) : {ct} {cmd.split()[1]} {cmd.split()[2]} {cmd.split()[6]} ...",end="")
                    system_command(cfg, cmd)

                    output_found = np.array([0, 0, 0])
                    for i, ext in enumerate(output_extensions):
                        outfn = f"{cfg.datevshot}_{cmd.split()[6]}{ext}"
                        if os.path.exists(os.path.join("detect_out/", outfn)):
                            output_found[i] = 1

                    if np.count_nonzero(output_found) != 3:
                        # something failed, we will want to re-run these once
                        #failed_list.append(cmd)
                        print(f" FAIL. Second attempt. No more retries.")
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

        print(f"continuum detections (rgetmax) ...")

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

        if os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.gaia"):
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


def amp_stats(cfg,shot_h5_fqfn=None):
    """

    this needs the shot h5 file to have been completed to work
    :param cfg:
    :return:
    """

    try:
        rc = 0
        if shot_h5_fqfn is None:
            shot_h5_fqfn = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")

        print(f"[{cfg.datevshot}] Computing amp statistics from: {shot_h5_fqfn} ... ")
        shot_dict = AmpStats.make_stats_for_shot(fqfn=shot_h5_fqfn,save=True,preload=False)


        if shot_dict is not None:
            t = AmpStats.stats_shot_dict_to_table(shot_dict)
            t = t[t['n_lo'] >= 0] #use n_lo column to select ... the -1 values are where this failed
                                  # (e.g. usually for dithers that don't exist)

            #???how much of stats_qc needs to be re-done since it is based on 3-dithers and some joint statistics???
            #several of the checks are looking for extreme variation over the dithers, which can't be done with just one dither
            t = AmpStats.stats_qc(t, extend=True)

            t.write(f"{cfg.datevshot}_ampstats.fits",format="fits",overwrite=True)
            t.write(f"{cfg.datevshot}_ampstats.tab", format="ascii",overwrite=True)

            #always creat the bad amps file, even if none trigger
            with open(f"{cfg.datevshot}_badamps.txt","w") as f:
                for row in t[t['flag']!=1]: #1 == Good, 0 = Bad
                    f.write(f"{row['multiframe']} exp{str(row['expnum']).zfill(2)} \n")


            #assuming these are post-HETDEX, go ahead and put this in the shot.h5 file
            #NOTICE: this is not done for the original HETDEX shots
            #needs the actual h5 file

            #NOTICE: the "flag" key DOES NOT EXIST here ... it is added to table t above, but not to the shot_dict
            print(f"Adding AmpStats table to  {shot_h5_fqfn} ...")
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
            rc = -1


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could produce amp statistics.", traceback.format_exc())

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

def shot_analyisis(cfg):
    """

    this needs the shot h5 file to have been completed to work
    :param cfg:
    :return:
    """

    try:
        print(f"[{cfg.datevshot}] Making basic IFU analysis images ...")
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

        analysis_dir = os.path.join(cfg.cwd,"analysis/ifus/")
        Path(analysis_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(analysis_dir)
        mfs_in_shot = np.unique(h5.root.Data.FiberIndex.read(field="multiframe"))
        ifus_in_shot = np.unique([mfs[0:-3] for mfs in mfs_in_shot])
        slotids = [ifu[10:13] for ifu in ifus_in_shot]
        #put these in order of slotid
        slotids, ifus_in_shot = zip(*sorted(zip(slotids, ifus_in_shot)))

        img = "clean_image"
        cmap = "gray"

        #just assume 3 dithers ... if they do not exist, they will be blank
        #there are a few that have 4 or more exposures and we will just ignore that
        for ii, mf_base in enumerate(ifus_in_shot):
            print(f"[{cfg.datevshot}] Making basic IFU analysis images: {mf_base.decode()}")
            plt.close('all')
            fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(9, 12))
            plot_config = list(np.arange(431, 443, 1))
            fig.suptitle(f"{cfg.datevshot} {mf_base.decode()}")

            for ai, amp in enumerate([b'_RU', b'_RL', b'_LL', b'_LU']):
                ei = 0  # exposure index
                # print(ai,ei,amp)
                mf = mf_base + amp
                data = h5.root.Data.Images.read_where("multiframe==mf", field=img)
                exp = h5.root.Data.Images.read_where("multiframe==mf", field='expnum')

                # ax = plt.subplot(plot_config[ci])
                ax = axes[ai, ei]
                ei += 1
                sel = exp == 1
                try:
                    if np.count_nonzero(sel) == 1:
                        ax.set_title(f"{amp.decode()} x1")
                        vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                        if amp != b'_LU':
                            ax.set_xticks([])
                        ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax)
                    else:
                        ax.set_xticks([])
                        ax.set_yticks([])

                except:
                    ax.set_xticks([])
                    ax.set_yticks([])

                sel = exp == 2
                ax = axes[ai, ei]
                ei += 1
                try:
                    if np.count_nonzero(sel) == 1:

                        ax.set_title(f"{amp.decode()} x2")
                        vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                        # ax.yaxis.label.set_visible(False)
                        ax.set_yticks([])
                        if amp != b'_LU':
                            ax.set_xticks([])
                        ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax)
                    else:
                        ax.set_xticks([])
                        ax.set_yticks([])

                except:
                    ax.set_xticks([])
                    ax.set_yticks([])

                sel = exp == 3
                ax = axes[ai, ei]
                ei += 1
                try:
                    if np.count_nonzero(sel) == 1:
                        ax.set_title(f"{amp.decode()} x3")
                        vmin, vmax = Utils.get_vrange(data[sel][0], contrast=0.25)
                        if amp != b'_LU':
                            ax.set_xticks([])
                        ax.set_yticks([])
                        ax.imshow(data[sel][0], cmap=cmap, vmin=vmin, vmax=vmax)
                    else:
                        ax.set_xticks([])
                        ax.set_yticks([])

                except:
                    ax.set_xticks([])
                    ax.set_yticks([])

            plt.tight_layout()
            plt.savefig(f"i{mf_base.decode()[10:13]}_{cfg.datevshot}_{mf_base.decode()}.png", dpi=96)


    except:
        #print(traceback.format_exc())
        rc = -1
        print(f"[{cfg.datevshot}] Could complete basic shot analysis.",traceback.format_exc())

    try:
        h5.close()
    except:
        pass

    return rc


# defunct now
# def build_detection_tables(cfg):
#     """
#
#     build a catalog astropy table (fits format) for the line and continuum detections
#
#     :param cfg:
#     :return:
#     """
#
#
#     #table definitions
#
#     #since there should never be more than maybe a few thousand lines at MOST (these are single observations)
#     #I don't think there is any need to break up the table creation into chunks and then vstack the chunks
#     #Lines Detections Table
#
#     mc_colnames = ['wave', 'wave_err', 'flux', 'flux_err', 'linewidth', 'linewidth_err',
#                    'continuum', 'continuum_err', 'sn', 'sn_err', 'chi2', 'chi2_err', 'ra', 'dec',
#                    'datevshot', 'noise_ratio', 'linewidth_fix', 'chi2_fix', 'chi2fib',
#                    'src_index', 'multiname', 'exp', 'xifu', 'yifu', 'xraw', 'yraw', 'weight',
#                    'apcor', 'sn_cen', 'flux_noise_1sigma', 'sn_3fib', 'sn_3fib_cen', 'dummy']
#
#     spec_colnames = ["wave1d", "spec1d_nc", "spec1d_nc_err", "counts1d", "counts1d_err",
#                      "apsum_counts", "apsum_counts_err", "dummy", "apcor", "flag_pix", "src_index"]
#
#     list_colnames = ["ra", "dec", "x_ifu", "y_ifu", "multiname", "expnum", "distance", "wave", "timestamp", "date",
#                      "obsid", "x_raw", "y_raw", "weight", "flag", "src_index"]
#
#     LT = Table(dtype=[
#         ('line_detectid', np.int64), #just a rolling ID number with a prefix
#         ('shotid', np.int64),        #YYYYMMDDSSS these should all be the same shotid, but might be handy if vstacking outside of this use
#         ('ifu','S11'),               #w/o leading "multi_"  so just "123_456_789" (e.g. the IFU)
#         ('amp','S2'),                # "LL,LU,RL,RU"
#         ('fibernum', 'S3'),          #"001"-"112"
#         ('expnum','S2'),             #"01" to "03" (usually)
#         ('ra', np.float32), ('dec', np.float32), #decimal degrees
#         ('xifu', np.float32), ('yifu', np.float32),  # decimal degrees
#         ('xraw', np.float32), ('yraw', np.float32),
#
#         ('sn', np.float32), ('sn_err', np.float32),
#         ('sn_cen', np.float32), ('sn_3fib', np.float32), ('sn_3fib_cen', np.float32),
#         ('noise_ratio', np.float32),
#         ('chi2', np.float32),  ('chi2_err', np.float32),  # of the fit
#         ('chi2fib', np.float32), ('chi2_fix', np.float32),
#         ('apcor', np.float32),
#         ('wave', np.float32), ('wave_err', np.float32),  # AA
#         ('linewidth', np.float32), ('linewidth_err', np.float32),  #fit sigma AA
#         ('lineflux',np.float32), ('lineflux_err',np.float32),  #of ths fit, in ergs/s/cm2
#         ('continuum', np.float32), ('continuum_err', np.float32),  # of ths fit, in ergs/s/cm2/AA
#
#
#         ('obs_fluxd', (np.float32, 1036)),      #local sky subtracted, 1d flux in 1e-17 ergs/s/cm2/AA (so /2AA) NOT dust corrected
#         ('obs_fluxd_err', (np.float32, 1036)),
#        # ('dust_x', (np.float32, 1036)),         #multiplier for dust correction
#        ])
#
#     #continuum detections table (few different columns)
#     CT = Table(dtype=[
#         ('cont_detectid', np.int64),  # just a rolling ID number with a prefix
#         ('shotid', np.int64),
#         # YYYYMMDDSSS these should all be the same shotid, but might be handy if vstacking outside of this use
#         ('ifu', 'S11'),  # w/o leading "multi_"  so just "123_456_789" (e.g. the IFU)
#         ('amp', 'S2'),  # "LL,LU,RL,RU"
#         ('fibernum', 'S3'),  # "001"-"112"
#         ('expnum', 'S2'),  # "01" to "03" (usually)
#         ('ra', np.float32), ('dec', np.float32),  # decimal degrees
#         ('xifu', np.float32), ('yifu', np.float32),  # decimal degrees
#         ('xraw', np.float32), ('yraw', np.float32),
#         ('apcor', np.float32),
#
#
#         ('obs_fluxd', (np.float32, 1036)), # local sky subtracted, 1d flux in 1e-17 ergs/s/cm2/AA (so /2AA) NOT dust corrected
#         ('obs_fluxd_err', (np.float32, 1036)),
#     ])
#
#     #get the bad amps
#
#     try:
#         #we are in the shot working dir
#         h5 = tables.open_file(f"{cfg.datevshot}.h5",mode='r')
#         bad_amps_list = list(h5.root.AmpStats.read_where("flag==0", field="multiframe").astype(str))
#         print(f"Loaded {len(bad_amps_list)} bad amps ...")
#         h5.close()
#     except:
#         print("Could not load bad_amps_list")
#         bad_amps_list = []
#
#     try:
#         print("Building lines catalog ...")
#         os.chdir(os.path.join(cfg.cwd))
#         mc_files = glob.glob("alldet/detect_out/*.mc")  # should be 1:1 with *.list and *.mc
#         datevshot = datevshot = cfg.datevshot
#         #todo: convert flux to fluxd (and flux_err to fluxd_err)
#         # e.g.    like         rowspectra["spec1d"] = dataspec["spec1d_nc"] / dataspec["apcor"]
#         #                      rowspectra["spec1d_err"] = dataspec["spec1d_nc_err"] / dataspec["apcor"]
#         #todo: apply apcor (divide fluxd and fluxd_err by apcor)
#         #todo: fetch and apply dust correction
#         #todo: fetch and apply OTHER corrections (like wd, ... see hetdex_api:create_detecet_hdf5.py
#
#         #todo: check on building shot h5 file ... see Google docs (ingestion) and hetdex_api
#
#         #Table is LT
#         T = LT
#         detectid_ct = np.int64(cfg.datevshot.replace('v','',))*np.int64(1e7) + 1000000  #YYYYMMDDSSS1000000
#
#         for mc in mc_files:
#
#             sp = mc.replace(".mc", ".spec")
#             ld = mc.replace(".mc", ".list")
#
#             if not os.path.exists(sp) or not os.path.exists(ld):
#                 continue
#
#             bn = os.path.basename(mc)
#             datevshot = str(bn[0:12])
#             shotid = np.int64(datevshot.replace('v', ''))
#             ifu = str(bn[13:24])  # spec_slot_ifuid"
#
#             try:
#                 t_mc = Table.read(mc, format="ascii.no_header", names=mc_colnames)
#                 t_sp = Table.read(sp, format="ascii.no_header", names=spec_colnames)
#                # t_ld = Table.read(ld, format="ascii.no_header", names=list_colnames)
#             except:
#                 continue
#
#
#             #multiple entries in each list (one per fiber) for each mc entry
#             for row in t_mc:
#                 src_index = row['src_index']
#
#                 t1 = t_mc[t_mc['src_index'] == src_index]  #one row
#                 t2 = t_sp[t_sp['src_index'] == src_index]  #many rows (one for each wavelength for the src_index)
#                 #don't actually need the *.list file for this version of the table
#                 #t3 = t_ld[t_ld['src_index'] == src_index]  #several rows, one for each fiber
#
#                 #sort t3 by weight so highest weight is top
#                 #t3.sort('weight').reverse()
#
#                 #check t1['multiname][0] vs bad amps
#                 if t1['multiname'][0][0:-4] in bad_amps_list:
#                     print(f"Skipping line detect in mc {mc} at index {src_index} for matching bad amp.")
#                     continue
#                 else: #for debug
#                     print(f"DEBUG: {t1['multiname'][0]} or {t1['multiname'][0][0:-4]} not in bad amps list.")
#
#                 detectid_ct += 1
#
#                 T.add_row([
#                     detectid_ct,
#                     shotid,
#                     ifu,
#                     t1['multiname'][0].split("_")[4], #amp #multi_323_043_040_LL_093
#                     t1['multiname'][0].split("_")[5],  # fibernum #multi_323_043_040_LL_093 (top fiber, should match sorted t3[0]
#                     #should be same as t3['multiname'][0] striped to fiber
#                     t1['exp'][0][-2:], #input is "exp01" just want the "01"
#                     t1['ra'][0],
#                     t1['dec'][0],
#                     t1['xifu'][0],
#                     t1['yifu'][0],
#                     t1['xraw'][0],
#                     t1['yraw'][0],
#                     t1['sn'][0],
#                     t1['sn_err'][0],
#                     t1['sn_cen'][0],
#                     t1['sn_3fib'][0],
#                     t1['sn_3fib_cen'][0],
#                     t1['noise_ratio'][0],
#                     t1['chi2'][0],
#                     t1['chi2_err'][0],
#                     t1['chi2fib'][0],
#                     t1['chi2_fix'][0],
#                     t1['apcor'][0],
#                     t1['wave'][0],
#                     t1['wave_err'][0],
#                     t1['linewidth'][0],
#                     t1['linewidth_err'][0],
#                     t1['flux'][0],
#                     t1['flux_err'][0],
#                     t1['continuum'][0],
#                     t1['continuum_err'][0],
#
#                     np.array(t2['spec1d_nc'] / 2.), # observed flux density (hence the /2.0 AA)
#                     np.array(t2['spec1d_nc_err'] / 2.),# error
#                 ])
#
#         tname = f"{datevshot}_line.fits"
#         T.write(tname,format="fits",overwrite=True)
#         print(f"Wrote raw lines table: {os.getcwd()}/{tname}")
#
#     except:
#         print(traceback.format_exc())
#         rc = -1
#         print(f"Exception building line detections table:  {cfg.datevshot}")
#
#     try:
#         print("Building continuum catalog ...")
#         os.chdir(os.path.join(cfg.cwd))
#         spec_files = glob.glob("cs/spec/*.spec") #should be 1:1 with *.list (there are no *.mc)
#
#         detectid_ct = np.int64(cfg.datevshot.replace('v', '', )) * np.int64(1e7) + 9000000
#         datevshot = cfg.datevshot
#
#         T = CT
#
#         for sp in spec_files:
#             #one each for each continuum detection (note: these are not clustered into sources)
#             ld = sp.replace(".spec", ".list")
#
#             if not os.path.exists(ld):
#                 continue
#
#             bn = os.path.basename(sp)
#             datevshot = str(bn[0:12])
#             shotid = np.int64(datevshot.replace('v', ''))
#             #ifu = str(bn[13:24])  # spec_slot_ifuid"
#
#             try:
#                 t_sp = Table.read(sp, format="ascii.no_header", names=spec_colnames) #one row for each wavelength bin
#                 t_ld = Table.read(ld, format="ascii.no_header", names=list_colnames) #one for for each fiber
#             except:
#                 continue
#
#             src_index = t_ld['src_index'][0] #all the same src_index ... usually this is a value of 1, but can be 0
#
#             # there is no *.mc for continuum sources
#             #t1 = t_mc[t_mc['src_index'] == src_index]  # one row
#             t2 = t_sp[t_sp['src_index'] == src_index]  # many rows (one for each wavelength for the src_index)
#             t3 = t_ld[t_ld['src_index'] == src_index]  #several rows, one for each fiber
#
#             # sort t3 by weight so the highest weight is top
#             # t3.sort('weight').reverse()
#
#
#             if t3['multiname'][0][0:-4] in bad_amps_list:
#                 print(f"Skipping continuum detect in sp {sp} at index {src_index} for matching bad amp.")
#                 continue
#
#             detectid_ct += 1
#
#             T.add_row([
#                     detectid_ct,
#                     shotid,
#                     t3['multiname'][0][6:17], #ifu
#                     t3['multiname'][0].split("_")[4], #amp #multi_323_043_040_LL_093
#                     t3['multiname'][0].split("_")[5][0:3],  # fibernum #multi_414_038_035_RL_083.ixy
#                     t3['expnum'][0][-2:], #input is "exp01" just want the "01"
#                     t3['ra'][0],
#                     t3['dec'][0],
#                     t3['x_ifu'][0],
#                     t3['y_ifu'][0],
#                     t3['x_raw'][0],
#                     t3['y_raw'][0],
#                     np.mean(t2['apcor']),
#
#
#                     np.array(t2['spec1d_nc'] / 2.), # observed flux density (hence the /2.0 AA)
#                     np.array(t2['spec1d_nc_err'] / 2.),# error
#                 ])
#
#         tname = f"{datevshot}_cont.fits"
#         T.write(tname, format="fits", overwrite=True)
#         print(f"Wrote raw continuum table: {os.getcwd()}/{tname}")
#
#     except:
#         print(traceback.format_exc())
#         rc = -1
#         print(f"Exception building line detections table:  {cfg.datevshot}")
# #end build_detection_tables


def check_detid_counter(cfg,cont,nominal_max):
    """

    :param cfg:
    :param cont: if True, check the continuum h5, else check the line h5 file
    :return:
    """

    try:
        if cont:
            which = "cont"
        else:
            which = "line"

        h5 = tables.open_file(f"{cfg.datevshot}_{which}.h5")
        ids = h5.root.Detections.read(field="detectid")
        maxid = np.max(ids)

        if maxid > nominal_max:
            rc = len(ids)
        else:
            rc = 0

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

        if len(line_tab) == 0: #there are no entries
            del line_tab
            print(f"[{cfg.datevshot}] WARNING! No line sources recored. Moving on to continuum.")
            rc = 1 #not fatal, but not a success
        else:
            line_h5 = tables.open_file(os.path.join(cfg.cwd, f"{cfg.datevshot}_line.h5"))

            #reminder: if fof clustering re-runs (the step before this one), the gmag will be wiped out
            if 'gmag' not in line_tab.columns:
                print(f"[{cfg.datevshot}] Computing gmags for Diagnose...")
                # subselect ... these do not currently have a gmag, so need to make one (since ELiXer has not run yet)
                line_tab['gmag'] = 99.0

                for i in range(len(line_tab)):
                    d = line_tab[i]['detectid']

                    rows = line_h5.root.Spectra.read_where("detectid==d")
                    if len(rows) != 1:
                        print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                        continue
                    row=rows[0]
                    gmag, gmag_unc, *_ = SU.get_best_gmag(row['spec1d']*1e-17,row['spec1d_err']*1e-17,G.CALFIB_WAVEGRID)
                    line_tab['gmag'][i] = gmag


                #update
                line_tab.write(name, format="ascii", overwrite=True)
                print(f"[{cfg.datevshot}] Updated lines table with gmag: {os.getcwd()}/{name}")

            # now select on gmag < 23
            sel = np.array(line_tab['gmag'] > 0.0) & np.array(line_tab['gmag'] < 23.0)
            line_tab = line_tab[sel]


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

            for i in range(len(line_tab)):
                d = line_tab[i]['detectid']

                rows = line_h5.root.Spectra.read_where("detectid==d")
                if len(rows) != 1:
                    print(f"Diagnose preselection spectra failure for {d}")
                    spec_2D.append(np.zeros(1036))
                    error_2D.append(np.zeros(1036))
                    apcor.append(np.zeros(1036))
                    continue

                row = rows[0]
                #line_tab['spec'][i] = rows['spec1d']
                #line_tab['spec_err'][i] = rows['spec1d_err']
                spec_2D.append(row['spec1d'])
                error_2D.append(row['spec1d_err'])
                apcor.append(row['apcor'])

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

            if 'gmag' not in cont_tab.columns:
                print(f"[{cfg.datevshot}] Computing gmags for Diagnose ...")
                # subselect ... these do not currently have a gmag, so need to make one (since ELiXer has not run yet)
                cont_tab['gmag'] = 99.0

                for i in range(len(cont_tab)):
                    d = cont_tab[i]['detectid']

                    rows = cont_h5.root.Spectra.read_where("detectid==d")
                    if len(rows) != 1:
                        print(f"[{cfg.datevshot}] Diagnose preselection gmag failure for {d}")
                        continue
                    row=rows[0]
                    gmag, gmag_unc, *_ = SU.get_best_gmag(row['spec1d']*1e-17,row['spec1d_err']*1e-17,G.CALFIB_WAVEGRID)
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

            for i in range(len(cont_tab)):
                d = cont_tab[i]['detectid']

                rows = cont_h5.root.Spectra.read_where("detectid==d")
                if len(rows) != 1:
                    print(f"[{cfg.datevshot}] Diagnose preselection spectra failure for {d}")
                    spec_2D.append(np.zeros(1036))
                    error_2D.append(np.zeros(1036))
                    apcor.append(np.zeros(1036))
                    continue

                row = rows[0]

                spec_2D.append(row['spec1d'])
                error_2D.append(row['spec1d_err'])
                apcor.append(row['apcor'])

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
        Path(elixdir).mkdir(parents=True, exist_ok=True)

        if FilterDetsOnBadAmps:
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
            bad_amps_list = []


        make_lines = True
        tab = Table.read(f"{cfg.datevshot}_line_sourcecat.tab",format="ascii")
        if len(tab) == 0:
            print(f"[{cfg.datevshot}] Error! no line sources recorded. ")
            make_lines = False
        else:
            sel = np.array([x['multiframe'] not in bad_amps_list for x in tab])
            print(f"[{cfg.datevshot}] Excluding {len(sel)-np.count_nonzero(sel)} / {len(sel)} line detections as residing on bad amps.")
            line_dets = list(tab['detectid'][sel])
            np.savetxt(os.path.join(elixdir, "line.dets"), line_dets, fmt="%d")
            del tab

        make_conts = True
        tab = Table.read(f"{cfg.datevshot}_cont_sourcecat.tab",format="ascii")
        if len(tab) == 0:
            print(f"[{cfg.datevshot}] Error! no continuum sources recorded. ")
            make_conts = False
        else:
            sel = np.array([x['multiframe'] not in bad_amps_list for x in tab])
            print(f"[{cfg.datevshot}] Excluding {len(sel) - np.count_nonzero(sel)} / {len(sel)} continuum detections as residing on bad amps.")
            cont_dets = list(tab['detectid'][sel])
            np.savetxt(os.path.join(elixdir, "cont.dets"), cont_dets, fmt="%d")
            del tab

        if make_lines or make_conts:
            shot_h5 = os.path.join(cfg.cwd,f"{cfg.datevshot}.h5")
            line_h5 = os.path.join(cfg.cwd,f"{cfg.datevshot}_line.h5")
            cont_h5 = os.path.join(cfg.cwd, f"{cfg.datevshot}_cont.h5")
            diagnose_tab = os.path.join(cfg.cwd,f"diagnose_classifications.tab")

            which_elixer = f"python {elixer_path}/selixer.py" #"selixer.test "
            elixer_base_cmd = f" -f --slurm 0 --nodes 1 --log info --shot_h5 {shot_h5} --diagnose {diagnose_tab} " \
                              f" --png --error 3.0 --neighborhood 10.0 --post_merge 2 "
            elixer_line_cmd = f" --name line --dets line.dets  --hdf5 {line_h5} "
            elixer_cont_cmd = f" --name cont --continuum --dets cont.dets  --hdf5 {cont_h5} "

            with open(os.path.join(elixdir,"elixer_cmd.txt"),"w") as f:
                if make_lines:
                    f.write(f"{which_elixer} {elixer_base_cmd} {elixer_line_cmd} \n")
                if make_conts:
                    f.write(f"{which_elixer} {elixer_base_cmd} {elixer_cont_cmd} \n")

            print(f"[{cfg.datevshot}]Preparing ELiXer SLURMS ... ")
            os.chdir(elixdir)
            system_command(cfg,"source elixer_cmd.txt")

            print("Default ELiXer SLURMS prepared. However, you may edit ./elixer/elixer_cmd.txt re-run if needed.")
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

if queue_elixer:
    run_queue_elixer(cfg)
    exit(0)

if cfg.clean_only:
    print(f"[{cfg.datevshot}] Performing only the CLEAN, level : {cfg.clean} ...")
    post_clean(cfg)
    Quit(cfg,0,"Clean complete. Exiting",write_status=False)

rc = precheck(cfg)
if rc < 0:
    Quit(cfg,rc,"Precheck failed. Reduction cannot run.",write_status=False)

cfg.numexp, cfg.gettar_fn = num_exposures_in_shot(cfg.shotid)

if cfg.numexp <= 0:
    Quit(cfg, -1, f"Could not find shot {cfg.datevshot}",write_status=False)


# if cfg.simul == 1:
#     NumProcs_mp_rcal = 10
#     NumProcs_mp_rf1 = 10

if cfg.exp <= 0:
    print(f"Working on {cfg.datevshot} with {cfg.numexp} exposure(s) ...")
    if cfg.numexp == 1 or cfg.numexp == 3: #okay
        pass
    else:
        #print(f"[{cfg.datevshot}] !!! bad !!! Unusual number of exposures ({cfg.numexp}).")
        print(f"********************************************************************************************")
        print(f"!!! WARNING !!! Unusual number of exposures ({cfg.numexp}) !!! Reduction may be problematic.")
        print(f"                You may want to consider reducing each exposure individually. ")
        print(f"********************************************************************************************")
else:
    if cfg.exp <= cfg.numexp:
        print(f"Working on {cfg.datevshot} exposure #{cfg.exp} ...")
        if int(cfg.datevshot[0:8]) < 20240800 and cfg.numexp == 3:
            try:
                dex_dvs = np.loadtxt("/corral-repl/utexas/Hobby-Eberly-Telesco/detect/fwhm.all",dtype=str,usecols=0)
                if cfg.datevshot in dex_dvs:
                    print(f"WARNING! {cfg.datevshot} is an original HETDEX observation.")
                    print(f"         Attempting to use only one exposure may not generate the expected results.")
                #else: we are not going to print anything ... while this occurred during the original HETDEX run
                #      this is not a previously reduced HETDEX observation
            except: #just assume it could be an original observation and warn
                print(f"WARNING! {cfg.datevshot} appears to be an original HETDEX observation.")
                print(f"         Attempting to use only one exposure may not generate the expected results.")
    else:
        Quit(cfg, -1, f"Invalid exposure. Requesting exp #{cfg.exp} but {cfg.datevshot} has only {cfg.numexp}",write_status=False)

rc = initial_setup(cfg)

if rc < 0:
    Quit(cfg,rc,"Could not complete initial setup.",write_status=False)

rc = node_setup(cfg)

if cfg.numexp < 3:
    print(f"Fewer than 3 exposures (assume dithers). Checking guider for seeing FWHM...")
    cfg.guider_fwhm = get_guider_fwhm(cfg)
    if cfg.guider_fwhm is not None:
        print(f"Using guider FWHM = {cfg.guider_fwhm}")
    else:
        print(f"Unable to obtain guider seeing FWHM. Will measure as best can be from available data.")




#########
# after the initial setup, move stdout and stderr to a log file
#########


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

###########
# step1
###########
if s01_run1s and not dtprog["s01_run1s"]:
    run_run1s(cfg)

    # todo: run any checks
    check_run1s(cfg)

    #todo: this would be manual here, I think, but CAN copy /red1/xxx to /scratch/local/projects
    #  all the various CoFe*.fits and multi*.fits ... these are also in the
    #  local d<shot><exp> folder in the two tar files (_co.tar and _mu.tar for the CoFe*.fits and multi*.fits respectively)

    progress_update(cfg,dtprog,"s01_run1s")
else:
    print(f"[{cfg.datevshot}] Skipping run1s")

###########
# step2
# VDRP
###########

if s02_vdrp and not dtprog["s02_vdrp"]:
    run_vdrp(cfg)

    check_vdrp(cfg)

    #todo: optional manual step here (need to be hetdex user), copy the *.dithall to
    #  /scratch/projects/hetdex/detect/dithall   (and /coral-repl/...)

    progress_update(cfg,dtprog, "s02_vdrp")

else:
    print(f"[{cfg.datevshot}] Skipping vdrp")



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
    if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.gaia")):
        star_cat_list.append("gaia")

    if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.sdss")):
        star_cat_list.append("sdss")

    if os.path.exists(os.path.join(cfg.cwd,f"vdrp/shifts/dithall.panstarrs")):
         star_cat_list.append("panstarrs")

    if len(star_cat_list) == 0:
        Quit(cfg, -1, "FATAL. Something wrong. No star catalogs available under vdrp/shifts.")


    for star_cat in star_cat_list:
        print(f"[{cfg.datevshot}] flux calibration: {star_cat}")
        run_fluxcalibration(cfg,star_cat)

        if check_fluxcalibration(cfg) < 0:
            print(f"[{cfg.datevshot}] flux calibration: {star_cat} failed. Trying next ... ")
        else:
            break #this one was good

    #todo: optional:  copy detect/tp/yyyymmddvssssedtp_f.dat  to /scratch/projects/hetdex/detectp and /corral-repl/xxx

    #todo: optional: update  /scratch/projects/hetdex/detect/fwhm.all and norm.all
    #                see update_fwhm_norm script
    #

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
    print(f"[{cfg.datevshot}] Skipping flux calibration")


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
        print(f"[{cfg.datevshot}] Skipping IFU centers ...")

    if s04b_rfft and not dtprog["s04b_rfft"]:
        print(f"[{cfg.datevshot}] Running rfft (this may take a while) ...")
        if run_rfft(cfg) != 0:
            Quit(cfg, -1, "FATAL. rfft fail. One or more expected outputs failed.")
        else:
            progress_update(cfg,dtprog, "s04b_rfft")
    else:
        print(f"[{cfg.datevshot}] Skipping rfft")

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
        print(f"[{cfg.datevshot}] Skipping rcal_all")

    if s04d_shot_h5 and not dtprog["s04d_shot_h5"]:
        rc = build_shot_h5(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. Could not build shot h5 file. Cannot continue with catalog creation.")
        progress_update(cfg,dtprog, "s04d_shot_h5")

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

    #basic shot analysis (mostly images for review)
    if s04f_analysis and not dtprog["s04f_analysis"]:
        rc = shot_analyisis(cfg)
        if rc < 0:
            print(f"[{cfg.datevshot}]Non-fatal. Could not complete basic shot analysis output.")
        progress_update(cfg,dtprog, "s04f_analysis")

    if dtprog["s04a_get_ifucens"] and dtprog["s04b_rfft"] and dtprog["s04c_rcal_all"] and dtprog["s04d_shot_h5"] and dtprog["s04e_amp_stats"]:
        progress_update(cfg,dtprog, "s04_make_shot")
else:
    print(f"[{cfg.datevshot}] Skipping sky subtraction")


###########
# step5
# line detections
# continuum detections
###########

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
        print("skipping rdet_rf1 (line detection)")

    if s05c_rgetmax  and not dtprog["s05c_rgetmax"]:
        print("Running rgetmax (continuum detection) ...")
        rc = rgetmax(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rgetmax fail.")
        elif rc > 0:
            print("rgetmax: Limited success. Non-fatal. Will continue")

        progress_update(cfg,dtprog, "s05c_rgetmax")
    else:
        print("skipping rgetmax (continuum detection)")

    if s05e_detection_hdf5 and not dtprog["s05e_detection_hdf5"]:
        rc = build_detection_hdf5(cfg)

    progress_update(cfg,dtprog, "s05e_detection_hdf5")

    if dtprog["s05b_rdet_rf1"] and dtprog["s05c_rgetmax"] and dtprog["s05e_detection_hdf5"]:
        progress_update(cfg,dtprog, "s05_detection")

else:
    print("Skipping detections")


##################################################
# step6
# combine detections
# run diagnose and elixer
# make source catalogs
#################################################

if s06_catalogs and not dtprog["s06_catalogs"]:

    print("Catalog creation ... ")

    if s06b_fof and not dtprog["s06b_fof"]:
        try:
            #line_tab = Table.read(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.fits"),format="fits")
            line_h5 = tables.open_file(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.h5"))
            line_tab = Table(line_h5.root.Detections.read())
            line_h5.close()

            #subselect "nominal" good
            esel = np.array(line_tab['continuum'] >= -3)
            esel = esel & np.array(line_tab['sn'] >= 4.8)
            esel = esel & np.array(line_tab['chi2'] <= 2.5)
            # this is a bit more liberal than standard HETDEX_API (1.6 and 14, I think)
            esel = esel & np.array(line_tab['linewidth'] >= 1.5) & np.array(line_tab['linewidth'] <= 16)
            esel = esel & np.array(line_tab['chi2fib'] <= 4.5)  # fairly restrictive, this is from .mc file column #19

            line_tab = line_tab[esel]

            if line_tab is not None:
                tname = f"{cfg.datevshot}_line_sourcecat.tab"
                line_tab.write(tname, format="ascii", overwrite=True)
                print(f"Initial lines source table: {os.getcwd()}/{tname}")

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
                        print(f"Updated lines source table: {os.getcwd()}/{tname}")

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
                print(f"[{cfg.datevshot}] Initial continuum source table: {os.getcwd()}/{tname}")

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
                        print(f"[{cfg.datevshot}] Updated continuum source table: {os.getcwd()}/{tname}")
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

    if s06d_elixer and not dtprog["s06d_elixer"]:

        rc = prep_elixer(cfg)

        progress_update(cfg,dtprog, "s06d_elixer")

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
