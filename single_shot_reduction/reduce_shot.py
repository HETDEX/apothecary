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
from dataclasses import dataclass

import tables
from astropy.table import Table, unique, vstack, join, Column, hstack, MaskedColumn
from h5tools import amp_stats as AmpStats
import hetdex_tools.fof_kdtree as fof


#just want the path for hetdex_api (see later)
import importlib.util

import traceback

########################################################################
# CONFIGURATION
########################################################################
EchoCmds = False #if True echo system commands to the log
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot
ElixerSnrThresh = 4.5 #do not run elixer on line sources where the S/N < 4.5

ScriptRepo = "/work/03261/polonius/hetdex/single_shot"
LocalScriptRepo = "./local_script_repo" #useful if running multiple single shots ... can copy remotely once
                                        #then copy locally from here for each shot
                                        #set to None if you do NOT want to use a local script dir cache
                                        #  and force a copy from the main repo each time

WorkDirRoot = "./"

HETRaw = "/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"
karlgettar = "/work/00115/gebhardt/maverick/gettar/"
karlfplane = "/work/00115/gebhardt/maverick/fplane/"
karlhome = "/home1/00115/gebhardt"
hetdex_api_path = os.path.dirname(importlib.util.find_spec("hetdex_api").origin)
#there is an extra "hetdex_api" at the end that points into the lower level directory for that. h5tools is actually a sibling
hetdex_api_path = "/".join(hetdex_api_path.split("/")[0:-1])


#execute steps
s01_run1s = True

s02_vdrp = True
do_panstarrs = False #only run PanSTARRS if true, otherwise just run the usual GAIA and SDSS

s03_fluxcal = True

s04_sky_subtraction = True
s04b_rfft = True
s04c_rcal_all = True
s04d_shot_h5 = True
s04e_amp_stats = True
s04_sky_subtraction = s04_sky_subtraction | s04b_rfft | s04c_rcal_all | s04d_shot_h5 | s04e_amp_stats #sanity catch

s05_detection = True
s05b_rdet_rf1 = True
s05c_rgetmax = True
s05d_detection_tables = True
s05e_detection_hdf5 = True
s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05d_detection_tables | s05e_detection_hdf5 #sanity catch


s06_catalogs = True
s06b_fof = True      #cluster the lines and continuum sources (separately)
s06c_diagnose = True #run Diagnose
s06d_elixer = True   #run elixer
s06e_source_cat = True #make a source catalog

s06_catalogs = s06_catalogs | s06b_fof | s06c_diagnose | s06d_elixer | s06e_source_cat

if False: #testing
    print("#################### TESTING ##########################")
    # execute steps
    s01_run1s = False

    s02_vdrp = False
    do_panstarrs = False  # only run PanSTARRS if true, otherwise just run the usual GAIA and SDSS

    s03_fluxcal = False

    s04_sky_subtraction = False
    s04b_rfft = False
    s04c_rcal_all = False
    s04d_shot_h5 = False
    s04e_amp_stats = False
    s04_sky_subtraction = s04_sky_subtraction | s04b_rfft | s04c_rcal_all | s04d_shot_h5 | s04e_amp_stats  # sanity catch

    s05_detection = False
    s05b_rdet_rf1 = False
    s05c_rgetmax = False
    s05d_detection_tables = False
    s05e_detection_hdf5 = False
    s05_detection = s05_detection | s05b_rdet_rf1 | s05c_rgetmax | s05d_detection_tables | s05e_detection_hdf5  # sanity catch

    s06_catalogs = True
    s06b_fof = True  # cluster the lines and continuum sources (separately)
    s06c_diagnose = False  # run Diagnose
    s06d_elixer = False  # run elixer
    s06e_source_cat = False  # make a source catalog

    s06_catalogs = s06_catalogs | s06b_fof | s06c_diagnose | s06d_elixer | s06e_source_cat





########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:

    clean: bool = False
    overwrite: bool = False
    resume: bool = False
    shotid: int = 0
    datevshot: str = ""
    exp: int = 0  #specific exposure number to reduce
    numexp: int = 0 #number of exposures in the shot
    cwd_orig: str = os.getcwd()
    cwd: str = os.getcwd()
    scriptdir: str = ""
    gettar_fn: str =  ""  # the runs* or runt* file from karlgettar folder with the date, shot, exp data
    starcat_ast = "gaia" #gaia, sdss, panstarrs   for Astrometry
    starcat_cal = "sdss" #gaia, sdss, panstarrs   for FluxCalibration

    orig_stdout = None
    orig_stderr = None
    file_stdout = None




########################################################################
# Basic user input
########################################################################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

cfg = Config()

if "-clean" in args:
    cfg.clean = True
    #idx = args.index("-clean")
    args.remove("-clean")

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
    try:
        #might have 'v' or 'd' or 's' as the separator between date and shot number
        cfg.shotid = int(args[0].replace("v","").replace("s","").replace("d",""))
    except:
        pass

    if not (20170101000 < cfg.shotid < FutureShotDateLimit):
        print(f"Fatal: Invalid shotid: {args[0]}")
        exit(-1)

    cfg.datevshot = str(cfg.shotid)[0:8] + "v" + str(cfg.shotid)[8:]



########################################################################
# worker functions
########################################################################

def Quit(cfg,rc,msg=None):
    """

    :param cfg:
    :param rc:
    :param msg:
    :return:
    """

    if msg is not None:
        print(f"({rc})",msg)
    else:
        print(f"({rc})")

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

    exit(rc)


def system_command(cfg,cmd):
    """

    wrapper to execute a system command

    :param cfg:
    :param cmd:
    :return:
    """

    #echo the command
    if EchoCmds:
        print("CMD: ", cmd)

    if cfg.file_stdout:
        os.system(f"{cmd} &>> {cfg.file_stdout.name}")
    else:
        os.system(f"{cmd}")


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
            print(f"Resuming. Leave directory intact: {workdir}")
            resume = True
        elif cfg.overwrite:
            print(f"Overwriting directory {workdir} ... ")
            shutil.rmtree(workdir)
        else:
            print(f"Shot directory already exists here! {workdir}")
            return -1

    if not resume:
        os.makedirs(workdir)

        if LocalScriptRepo is not None:
            if os.path.exists(LocalScriptRepo): #we want to use it
                print("Using existing local repo ...")
                cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
            else:
                #copy first to local script repo
                print("Copying to local repo ...")
                shutil.copytree(os.path.join(ScriptRepo, "science_reductions"),
                                os.path.join(os.getcwd(),LocalScriptRepo), dirs_exist_ok=True)
                cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
        else:
            print("Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")
    else:
        if LocalScriptRepo is not None:
            if os.path.exists(LocalScriptRepo): #we want to use it
                print("Using existing local repo ...")
                cfg.scriptdir = os.path.join(os.getcwd(), LocalScriptRepo)
            else:
                print("Fatal! --resume selected, but no script repo.")
                return -1
        else:
            print("Using main script repo (may be remote) ...")
            cfg.scriptdir = os.path.join(ScriptRepo,"science_reductions")


    os.chdir(workdir)
    cfg.cwd = os.getcwd() #now under the sci<shot> directory


    if not resume:
        print("Copying source code ...")
        ## if ANY of this fails it is fatal

        #shutil.copy2(os.path.join(cfg.scriptdir, "science_reductions", "rsetups"),".") #no, this function is its equivalent

        shutil.copy2(os.path.join(cfg.scriptdir, "rfixspec"), ".")
        shutil.copytree(os.path.join(cfg.scriptdir, "sciscripts"), ".", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir,"vdrp"), "vdrp", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "detect"), "detect", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "getcen"), "getcen", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "alldet"), "alldet", dirs_exist_ok=True)
        shutil.copytree(os.path.join(cfg.scriptdir, "cs"), "cs", dirs_exist_ok=True)

        #update the "home" path tilde
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits")
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_fix")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_field")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rback_fix")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbacks")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rbfits_s")  # use '#' as sed separator rather than "/"
        system_command(cfg,f"sed -i s#~gebhardt#{karlhome}# rimarb")  # use '#' as sed separator rather than "/"


        #update the old red1 path
        #yes, I want cwd_orig ... I want a single "reductions" direcetory off the top with all the sciXXX as siblings
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_field")  # use '#' as sed separator rather than "/"
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rback_fix")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rerun2")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# rtaremc") #not necessary, just for completeness
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtar")
        system_command(cfg, f"sed -i s#/scratch/03261/polonius/red1#{cfg.cwd_orig}# runtarm.defunct") #not necessary, just for completeness


        #extra files needed
        #vdrp : /work/00115/gebhardt/maverick/fplane ... need the fplane for the date
        os.makedirs(os.path.join(cfg.cwd,"vdrp/fplane"), exist_ok=True)
        shutil.copy2(os.path.join(karlfplane, f"fp{cfg.datevshot[0:8]}"), os.path.join(cfg.cwd,"vdrp/fplane"))

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
        rta_date, rta_shot, rta_ra, rta_dec, rta_v = np.loadtxt(os.path.join(karlgettar,f"rta.{cfg.datevshot[0:6]}"),
                                                                usecols=[1,2,3,4,5],unpack=True,dtype=str)
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

    print("check shout ifu ... ")
    try:
        wildcard = cfg.datevshot[0:6]
    except:
        wildcard = ""

    paths = glob.glob(f"{path}{wildcard}*v???")

    print(f"Checking {len(paths)} directories ... ")

    # for path in tqdm(paths):
    for path in paths:
        basedir = None
        try:
            fn = os.path.join(path, "shout.ifu")

            basedir = os.path.abspath(fn).split("/")[-2]

            stats = os.stat(fn)

            if stats.st_size == 0:
                print(f"{basedir} : empty shout.ifu")
                rc = -1
            elif stats.st_size < 1000:
                print(f"{basedir} : small shout.ifu ({stats.st_size})")
                rc = -1

        except Exception as E:
            print(f"{basedir} : unknown or missing shout.ifu")

    if rc != 0:
        print("check shout ifu ... fail ")
    else:
        print("check shout ifu ... OK")
    return rc

def vdrp_cp2dithall(cfg,catalog=None):
    """

    :param cfg:
    :return:
    """

    print("cp2dithall ... ")
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
                print(f"cp2dithall {catalog} ... fail. File does not exist {os.getcwd()} {dithall_use}")
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
        except Exception as E:
            print(E)
            rc = -1

    if rc != 0:
        print(f"cp2dithall {catalog}... fail")
    else:
        print(f"cp2dithall {catalog}... OK")
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

    print("VDRP: GAIA")
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
        print(f"VDRP: GAIA fail.",e,  "\n", traceback.format_exc())


    #!!! notice: we are doing limited checking here ... that will be done later
    # previously would run: check_norms YYYYMM   (included)
    #                       check_shot.ifu YYYYMM (included)
    #  examine the .pngs manually  NOT DONE
    #                  run: make_good_shots YYYYMM  NOT DONE
    #                       cp2dithall YYYYMM  (included)


    #SDSS
    print("VDRP: SDSS")
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
        print(f"VDRP: SDSS fail.", e, "\n", traceback.format_exc())


    if do_panstarrs:
        #PanSTARRS
        print("VDRP: PANSTARRS")
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
            print(f"VDRP: PANSTARRS fail.", e, "\n", traceback.format_exc())

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



def run_fluxcalibration(cfg,star_catalog='sdss'):
    """

    :param cfg:
    :return:
    """
    print("flux calibration (rallcal) ... ")

    #setup the softlink for the star_catalog
    os.chdir(os.path.join(cfg.cwd, "vdrp/shifts"))
    if os.path.exists(cfg.datevshot):
        system_command(cfg, f"unlink {cfg.datevshot}")
    if os.path.exists("dithall"):
        system_command(cfg, "unlink dithall")

    system_command(cfg, f"ln -s {star_catalog}/{cfg.datevshot} {cfg.datevshot}")
    system_command(cfg, f"ln -s dithall.{star_catalog} dithall")

    os.chdir(os.path.join(cfg.cwd,"detect"))


    #call rsetstar independently
    system_command(cfg, f"rsetstar {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} {star_catalog}")


    #no longer includes rsetstar
    system_command(cfg,f"rallcal {cfg.datevshot[0:8]} {cfg.datevshot[-3:]}")


    return 0


def check_fluxcalibration(cfg):
    """

    :param cfg:
    :return:
    """
    print("todo: check flux calibration ... ")
    rc = 0
    #todo: if tp/*setp_.dat is "bad" (last columns all zero)
    #      the try with a different star catalog ... e.g.
    #      by default, SDSS is used for calibration, but if that fails
    #      try GAIA and if that fails
    #      try PANSTARRS
    #          >> to do that, need to change the softlink under vdrp/shifts/YYYYMMDDsSSS to gaia/YYYYMMDDvSSS or panstarrs/YYYYMMDDvSSS
    #          >> if they exist (e.g. that might have failed during vdrp)
    #          >> then, rerun   run_fluxcalibration() and check again

    try:
        #20240730v006sedtp_f.dat
        out = np.loadtxt(f"tp/{cfg.datevshot}sedtp_f.dat")
        if not np.any(out[:, 5]):  # 5 is the actual throughput, I think and 4 is tied to it? cols 2, 3 don't seem to matter
            rc = -1
            print(f"bad throughput {cfg.datevshot}")

    except Exception as E:
        print(E)
        rc = -1
        print(f"bad throughput {cfg.datevshot}")

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
            #all is good, assume anyway
            rc = 0
        else:
            rc = -1

    except Exception as E:
        print(E)
        rc = -1
        print(f"failed to get IFU centers:  {cfg.datevshot}")


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
                    print(f"[Pass] {fn}")
                else:
                    print(f"[FAIL] {fn} file not found")
                    rc = -1 #even though this is fatal, go ahead and loop over all files and expXX so can get into the log
                            #could be a useful diagnoistic

    except Exception as E:
        print(E)
        rc = -1
        print(f"failed to get IFU centers:  {cfg.datevshot}")


    return rc



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

        multis = np.loadtxt(os.path.join("../getcen/",f"ifucen_{cfg.datevshot}.dat"),dtype=str,usecols=(0),unpack=True)
        ras,decs = np.loadtxt(os.path.join("../getcen/", f"ifucen_{cfg.datevshot}.dat"), dtype=float,usecols=(1,2),unpack=True)
        print(f"run_rcal ({len(ras)})...")
        ct = 0
        for multi,ra,dec in zip(multis,ras,decs):
            ct +=1
            print(f"{ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ",end="")
            system_command(cfg,f"rcal_all {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 106")

            #check the output exists cal_out/20240730v009_514_103_019_cal.fits
            outfn = f"{cfg.datevshot}_{multi[6:]}_cal.fits"
            if os.path.exists(os.path.join("cal_out/", outfn)):
                print("pass")
                passed_rcal_list.append(outfn)
            else:
                print("FAIL")
                failed_rcal_list.append(outfn)

        if len(passed_rcal_list) == len(ras):
            rc = 0
            print("All pass")
        elif len(failed_rcal_list) == len(ras):  # all failed
            rc = -1
            print("ALL failed")
        else:
            rc = 1
            print(f"Mixed results of {len(ras)}: {len(passed_rcal_list)} pass, {len(failed_rcal_list)} FAIL")


    except Exception as E:
        print(E)
        rc = -1
        print(f"Fatal exception in run_rcal:  {cfg.datevshot}")



    return rc


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

        print(f"line detections (rdet_rf1) ({len(ras)})...")
        ct = 0
        for multi, ra, dec in zip(multis, ras, decs):
            ct +=1
            print(f"{ct}) {multi[6:]}  {ra:0.7f} {dec:0.7f} ... ",end="")
            cmd = f"rf1 {ra:0.7f} {dec:0.7f} 35 4505 50 {multi[6:]} {cfg.datevshot} 1.70 3.0 3.5 0.5 3 104\n"
            system_command(cfg,cmd)

            # check the output exists cal_out/20240730v009_514_103_019_cal.fits

            #todo: assume good for now ... need to check
            #print("done (todo: check output files)")
            #check that .list, .mc, .spec exist
            #20240730v009_025_067_032.list

            output_found = np.array([0,0,0])
            for i,ext in enumerate(output_extensions):
                outfn = f"{cfg.datevshot}_{multi[6:]}{ext}"

                if os.path.exists(os.path.join("detect_out/", outfn)):
                    output_found[i] = 1

            if np.count_nonzero(output_found) != 3:
                #something failed, we will want to re-run these once
                failed_list.append(cmd)
                print(f"FAIL. May re-run at the end.")
            else:
                print(f"pass")


        if len(failed_list) > 0:
            if len(failed_list) < len(ras):
                print(f"{len(failed_list)} failed. Can be transient issues, so will re-run ...")
                for ct,cmd in enumerate(failed_list):
                    print(f"{ct} {cmd.split()[1]} {cmd.split()[2]} {cmd.split()[6]} ...",end="")
                    system_command(cfg, cmd)

                    output_found = np.array([0, 0, 0])
                    for i, ext in enumerate(output_extensions):

                        outfn = f"{cfg.datevshot}_{cmd.split()[6]}{ext}"

                        if os.path.exists(os.path.join("detect_out/", outfn)):
                            output_found[i] = 1

                    if np.count_nonzero(output_found) != 3:
                        # something failed, we will want to re-run these once
                        #failed_list.append(cmd)
                        print(f"FAIL. Second attempt. No more retries.")
                        rc = 1 #some failures, but not all
                    else:
                        print(f"pass")


            else:
                print(f"All failed. Will not attempt full re-run.")


    except Exception as E:
        print(E)
        rc = -1
        print(f"Fatal exception in rdet_rf1:  {cfg.datevshot}")

    return rc



#continuum detection
def rgetmax(cfg):
    """
    continuum detections

    :param cfg:
    :return:
    """

    try:


        failed_list = []
        #passed_list = []
        output_extensions = [".list", ".mc", ".spec"]

        rc = 0
        os.chdir(os.path.join(cfg.cwd, "cs"))  # make sure we are in the right directory
        os.makedirs("spec", exist_ok=True)

        print(f"continuum detections (rgetmax) ...")
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
            print(f"Continuum detections fail. Not fatal. cs/spec/{cfg.datevshot}cs.tar not created.")


    except Exception as E:
        print(E)
        rc = -1
        print(f"Fatal exception in rdet_rf1:  {cfg.datevshot}")

    return rc



def build_shot_h5(cfg):
    """

    :param cfg:
    :return:
    """

    try:
        print("Constructing shot hdf5 file. This will take a while ... ")
        rc = 0
        os.chdir(os.path.join(cfg.cwd))

        #needed for hetdex_api
        os.makedirs("match_pngs",exist_ok=True)

        #########################
        # initial hdf5 file
        ########################

        cmd = f"python3 {hetdex_api_path}/h5tools/create_shot_hdf5.py"
        cmd += " --tar"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}.h5\""
        cmd += f" --rootdir \"{cfg.cwd}\""

        if os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.gaia"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.gaia\""
        elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.sdss"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.sdss\""
        elif os.path.exists(f"{cfg.cwd}/vdrp/shifts/dithall.panstarrs"):
            cmd += f" --detect_path \"{cfg.cwd}/vdrp/shifts/dithall.panstarrs\""
        else:
            print("Fatal: cannot find *.dithall file for this shot")
            return -1



        system_command(cfg, cmd)

        #assume good?
        if not os.path.exists(f"{cfg.datevshot}.h5"):
            rc = -1
            return rc

        print(f"Created: {cfg.cwd}/{cfg.datevshot}.h5")

        #########################
        # now append_calfib
        ########################
        print("Appending calibrated fibers ... ")
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
        print("Appending fullsky model ... ")
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
        print("Create cal hdf5  ... ")


        #hetdex_api needs the local fwhm.all in a different format
        #detect/20240730v009
        if not os.path.exists(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail"):
            system_command(cfg,f"mv {cfg.cwd}/detect/{cfg.datevshot}/fwhm.all {cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail")
        fwhm, err, ns = np.loadtxt(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.detail",unpack=True,usecols=[0,1,2],max_rows=1,dtype=float)
        with open(f"{cfg.cwd}/detect/{cfg.datevshot}/fwhm.all","w") as f:
            f.write(f"{cfg.datevshot} {fwhm} {err} {int(ns)}\n")


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
        print("Create astrometry  ... ")
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
            print("Fatal: cannot find *.dithall file for this shot")
            return -1

        system_command(cfg, cmd)

        # #for convenience with legacy functions, duplicate the Shot group as a Survey group
        # shot_h5 = tables.open_file(f"{cfg.datevshot}.h5","a")
        # shot_h5.root.Shot._f_copy(newparent="root", newname="Survey", recursive=True, createparents=True)
        # shot_h5.close()

        print(f"Done: {cfg.cwd}/{cfg.datevshot}.h5")

    except Exception as E:
        print(E)
        rc = -1
        print(f"Could not build hdf5 shot file for  {cfg.datevshot}")

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

        print(f"Computing amp statistics from: {shot_h5_fqfn} ... ")
        shot_dict = AmpStats.make_stats_for_shot(fqfn=shot_h5_fqfn,save=True,preload=False)

        # if shot_dict is not None:
        #     #... not sure I want to actually modify the h5 file
        #     # the code is a bit dated (in hetdex api) and we had discussions about leaving the shot h5 files alone
        #     # ... this, though, is different now since these are not HETDEX and are working on individual shots
        #     #needs the actual h5 file
        #     h5 = tables.open_file(shot_h5_fqfn,mode="a")
        #     AmpStats.stats_update_shot(h5,shot_dict)
        # else:
        #     print(f"FAIL. Could not compute amp stats.")
        #     rc = -1


        if shot_dict is not None:
            t = AmpStats.stats_shot_dict_to_table(shot_dict)
            t = t[t['n_lo'] >= 0] #use n_lo column to select ... the -1 values are where this failed
                                  # (e.g. usually for dithers that don't exist)
            #???how much of stats_qc needs to be re-done since it is based on 3-dithers and some joint statistics???
            #several of the checks are looking for extreme variation over the dithers, which can't be done with just one dither
            t = AmpStats.stats_qc(t, extend=True)

            t.write(f"{cfg.datevshot}_ampstats.fits",format="fits")
            t.write(f"{cfg.datevshot}_ampstats.tab", format="ascii")

            #always creat the bad amps file, even if none trigger
            with open(f"{cfg.datevshot}_badamps.txt","w") as f:
                for row in t[t['flag']>0]:
                    f.write(f"{row['multiframe']} exp{str(row['expnum']).zfill(2)} \n")



        else:
            print(f"FAIL. Could not compute amp stats.")
            rc = -1


    except Exception as E:
        print(E)
        rc = -1
        print(f"Could produce amp statistics for:  {cfg.datevshot}")

    return rc

def build_detection_tables(cfg):
    """

    build a catalog astropy table (fits format) for the line and continuum detections

    :param cfg:
    :return:
    """


    #table definitions

    #since there should never be more than maybe a few thousand lines at MOST (these are single observations)
    #I don't think there is any need to break up the table creation into chunks and then vstack the chunks
    #Lines Detections Table

    mc_colnames = ['wave', 'wave_err', 'flux', 'flux_err', 'linewidth', 'linewidth_err',
                   'continuum', 'continuum_err', 'sn', 'sn_err', 'chi2', 'chi2_err', 'ra', 'dec',
                   'datevshot', 'noise_ratio', 'linewidth_fix', 'chi2_fix', 'chi2fib',
                   'src_index', 'multiname', 'exp', 'xifu', 'yifu', 'xraw', 'yraw', 'weight',
                   'apcor', 'sn_cen', 'flux_noise_1sigma', 'sn_3fib', 'sn_3fib_cen', 'dummy']

    spec_colnames = ["wave1d", "spec1d_nc", "spec1d_nc_err", "counts1d", "counts1d_err",
                     "apsum_counts", "apsum_counts_err", "dummy", "apcor", "flag_pix", "src_index"]

    list_colnames = ["ra", "dec", "x_ifu", "y_ifu", "multiname", "expnum", "distance", "wave", "timestamp", "date",
                     "obsid", "x_raw", "y_raw", "weight", "flag", "src_index"]

    LT = Table(dtype=[
        ('line_detectid', np.int64), #just a rolling ID number with a prefix
        ('shotid', np.int64),        #YYYYMMDDSSS these should all be the same shotid, but might be handy if vstacking outside of this use
        ('ifu','S11'),               #w/o leading "multi_"  so just "123_456_789" (e.g. the IFU)
        ('amp','S2'),                # "LL,LU,RL,RU"
        ('fibernum', 'S3'),          #"001"-"112"
        ('expnum','S2'),             #"01" to "03" (usually)
        ('ra', np.float32), ('dec', np.float32), #decimal degrees
        ('xifu', np.float32), ('yifu', np.float32),  # decimal degrees
        ('xraw', np.float32), ('yraw', np.float32),

        ('sn', np.float32), ('sn_err', np.float32),
        ('sn_cen', np.float32), ('sn_3fib', np.float32), ('sn_3fib_cen', np.float32),
        ('noise_ratio', np.float32),
        ('chi2', np.float32),  ('chi2_err', np.float32),  # of the fit
        ('chi2fib', np.float32), ('chi2_fix', np.float32),
        ('apcor', np.float32),
        ('wave', np.float32), ('wave_err', np.float32),  # AA
        ('linewidth', np.float32), ('linewidth_err', np.float32),  #fit sigma AA
        ('lineflux',np.float32), ('lineflux_err',np.float32),  #of ths fit, in ergs/s/cm2
        ('continuum', np.float32), ('continuum_err', np.float32),  # of ths fit, in ergs/s/cm2/AA


        ('obs_fluxd', (np.float32, 1036)),      #local sky subtracted, 1d flux in 1e-17 ergs/s/cm2/AA (so /2AA) NOT dust corrected
        ('obs_fluxd_err', (np.float32, 1036)),
       # ('dust_x', (np.float32, 1036)),         #multiplier for dust correction
       ])

    #continuum detections table (few different columns)
    CT = Table(dtype=[
        ('cont_detectid', np.int64),  # just a rolling ID number with a prefix
        ('shotid', np.int64),
        # YYYYMMDDSSS these should all be the same shotid, but might be handy if vstacking outside of this use
        ('ifu', 'S11'),  # w/o leading "multi_"  so just "123_456_789" (e.g. the IFU)
        ('amp', 'S2'),  # "LL,LU,RL,RU"
        ('fibernum', 'S3'),  # "001"-"112"
        ('expnum', 'S2'),  # "01" to "03" (usually)
        ('ra', np.float32), ('dec', np.float32),  # decimal degrees
        ('xifu', np.float32), ('yifu', np.float32),  # decimal degrees
        ('xraw', np.float32), ('yraw', np.float32),
        ('apcor', np.float32),


        ('obs_fluxd', (np.float32, 1036)), # local sky subtracted, 1d flux in 1e-17 ergs/s/cm2/AA (so /2AA) NOT dust corrected
        ('obs_fluxd_err', (np.float32, 1036)),
    ])

    try:
        print("Building lines catalog ...")
        os.chdir(os.path.join(cfg.cwd))
        mc_files = glob.glob("alldet/detect_out/*.mc")  # should be 1:1 with *.list and *.mc
        datevshot = datevshot = cfg.datevshot
        #todo: convert flux to fluxd (and flux_err to fluxd_err)
        # e.g.    like         rowspectra["spec1d"] = dataspec["spec1d_nc"] / dataspec["apcor"]
        #                      rowspectra["spec1d_err"] = dataspec["spec1d_nc_err"] / dataspec["apcor"]
        #todo: apply apcor (divide fluxd and fluxd_err by apcor)
        #todo: fetch and apply dust correction
        #todo: fetch and apply OTHER corrections (like wd, ... see hetdex_api:create_detecet_hdf5.py

        #todo: check on building shot h5 file ... see Google docs (ingestion) and hetdex_api

        #Table is LT
        T = LT
        detectid_ct = np.int64(cfg.datevshot.replace('v','',))*np.int64(1e7) + 1000000  #YYYYMMDDSSS1000000

        for mc in mc_files:

            sp = mc.replace(".mc", ".spec")
            ld = mc.replace(".mc", ".list")

            if not os.path.exists(sp) or not os.path.exists(ld):
                continue

            bn = os.path.basename(mc)
            datevshot = str(bn[0:12])
            shotid = np.int64(datevshot.replace('v', ''))
            ifu = str(bn[13:24])  # spec_slot_ifuid"

            try:
                t_mc = Table.read(mc, format="ascii.no_header", names=mc_colnames)
                t_sp = Table.read(sp, format="ascii.no_header", names=spec_colnames)
               # t_ld = Table.read(ld, format="ascii.no_header", names=list_colnames)
            except:
                continue


            #multiple entries in each list (one per fiber) for each mc entry
            for row in t_mc:
                src_index = row['src_index']

                t1 = t_mc[t_mc['src_index'] == src_index]  #one row
                t2 = t_sp[t_sp['src_index'] == src_index]  #many rows (one for each wavelength for the src_index)
                #don't actually need the *.list file for this version of the table
                #t3 = t_ld[t_ld['src_index'] == src_index]  #several rows, one for each fiber

                #sort t3 by weight so highest weight is top
                #t3.sort('weight').reverse()

                detectid_ct += 1

                T.add_row([
                    detectid_ct,
                    shotid,
                    ifu,
                    t1['multiname'][0].split("_")[4], #amp #multi_323_043_040_LL_093
                    t1['multiname'][0].split("_")[5],  # fibernum #multi_323_043_040_LL_093 (top fiber, should match sorted t3[0]
                    #should be same as t3['multiname'][0] striped to fiber
                    t1['exp'][0][-2:], #input is "exp01" just want the "01"
                    t1['ra'][0],
                    t1['dec'][0],
                    t1['xifu'][0],
                    t1['yifu'][0],
                    t1['xraw'][0],
                    t1['yraw'][0],
                    t1['sn'][0],
                    t1['sn_err'][0],
                    t1['sn_cen'][0],
                    t1['sn_3fib'][0],
                    t1['sn_3fib_cen'][0],
                    t1['noise_ratio'][0],
                    t1['chi2'][0],
                    t1['chi2_err'][0],
                    t1['chi2fib'][0],
                    t1['chi2_fix'][0],
                    t1['apcor'][0],
                    t1['wave'][0],
                    t1['wave_err'][0],
                    t1['linewidth'][0],
                    t1['linewidth_err'][0],
                    t1['flux'][0],
                    t1['flux_err'][0],
                    t1['continuum'][0],
                    t1['continuum_err'][0],

                    np.array(t2['spec1d_nc'] / 2.), # observed flux density (hence the /2.0 AA)
                    np.array(t2['spec1d_nc_err'] / 2.),# error
                ])

        tname = f"{datevshot}_line.fits"
        T.write(tname,format="fits",overwrite=True)
        print(f"Wrote raw lines table: {os.getcwd()}/{tname}")

    except Exception as E:
        print(E)
        rc = -1
        print(f"Exception building line detections table:  {cfg.datevshot}")

    try:
        print("Building continuum catalog ...")
        os.chdir(os.path.join(cfg.cwd))
        spec_files = glob.glob("cs/spec/*.spec") #should be 1:1 with *.list (there are no *.mc)

        detectid_ct = np.int64(cfg.datevshot.replace('v', '', )) * np.int64(1e7) + 9000000
        datevshot = cfg.datevshot

        T = CT

        for sp in spec_files:
            #one each for each continuum detection (note: these are not clustered into sources)
            ld = sp.replace(".spec", ".list")

            if not os.path.exists(ld):
                continue

            bn = os.path.basename(sp)
            datevshot = str(bn[0:12])
            shotid = np.int64(datevshot.replace('v', ''))
            #ifu = str(bn[13:24])  # spec_slot_ifuid"

            try:
                t_sp = Table.read(sp, format="ascii.no_header", names=spec_colnames) #one row for each wavelength bin
                t_ld = Table.read(ld, format="ascii.no_header", names=list_colnames) #one for for each fiber
            except:
                continue

            src_index = t_ld['src_index'][0] #all the same src_index ... usually this is a value of 1, but can be 0

            # there is no *.mc for continuum sources
            #t1 = t_mc[t_mc['src_index'] == src_index]  # one row
            t2 = t_sp[t_sp['src_index'] == src_index]  # many rows (one for each wavelength for the src_index)
            t3 = t_ld[t_ld['src_index'] == src_index]  #several rows, one for each fiber

            # sort t3 by weight so the highest weight is top
            # t3.sort('weight').reverse()

            detectid_ct += 1

            T.add_row([
                    detectid_ct,
                    shotid,
                    t3['multiname'][0][6:17], #ifu
                    t3['multiname'][0].split("_")[4], #amp #multi_323_043_040_LL_093
                    t3['multiname'][0].split("_")[5][0:3],  # fibernum #multi_414_038_035_RL_083.ixy
                    t3['expnum'][0][-2:], #input is "exp01" just want the "01"
                    t3['ra'][0],
                    t3['dec'][0],
                    t3['x_ifu'][0],
                    t3['y_ifu'][0],
                    t3['x_raw'][0],
                    t3['y_raw'][0],
                    np.mean(t2['apcor']),


                    np.array(t2['spec1d_nc'] / 2.), # observed flux density (hence the /2.0 AA)
                    np.array(t2['spec1d_nc_err'] / 2.),# error
                ])

        tname = f"{datevshot}_cont.fits"
        T.write(tname, format="fits", overwrite=True)
        print(f"Wrote raw continuum table: {os.getcwd()}/{tname}")

    except Exception as E:
        print(E)
        rc = -1
        print(f"Exception building line detections table:  {cfg.datevshot}")
#end build_detection_tables


def build_detection_hdf5(cfg):
    """

    line detects and continuum

    build a catalog astropy table (fits format) for the line and continuum detections

    :param cfg:
    :return:
    """


    try:

        detectid_base = np.int64(cfg.datevshot.replace('v', '', )) * np.int64(1e7) + 1000000
        print("Building detections hdf5 ... ")
        #cmd = f"python3 {hetdex_api_path}/h5tools/create_detect_hdf5.py"
        cmd = f"python3 {cfg.cwd}/../create_detect_hdf5.py"
        #cmd += f" --survey NA"
        cmd += f" --detectid_base {detectid_base}"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}_line.h5\""
        cmd += f" --detect_path \"{cfg.cwd}/alldet/detect_out\""

        system_command(cfg, cmd)

        detectid_base = np.int64(cfg.datevshot.replace('v', '', )) * np.int64(1e7) + 9000000
        #cmd = f"python3 {hetdex_api_path}/h5tools/create_cont_hdf5.py"
        cmd = f"python3 {cfg.cwd}/../create_cont_hdf5.py"
        cmd += f" --detectid_base {detectid_base}"
        cmd += f" --date {cfg.datevshot[0:8]}"
        cmd += f" --observation \"{cfg.datevshot[-3:]}\""
        cmd += f" -of \"{cfg.datevshot}_cont.h5\""
        cmd += f" --detect_path \"{cfg.cwd}/cs/spec\""

        system_command(cfg, cmd)

    except Exception as E:
        print(E)
        rc = -1
        print(f"Exception building line detections table:  {cfg.datevshot}")

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


cfg.numexp, cfg.gettar_fn = num_exposures_in_shot(cfg.shotid)

if cfg.numexp <= 0:
    Quit(cfg, -1, f"Could not find shot {cfg.datevshot}")

if cfg.exp <= 0:
    print(f"Working on {cfg.datevshot} with {cfg.numexp} exposure(s) ...")
else:
    if cfg.exp <= cfg.numexp:
        print(f"Working on {cfg.datevshot} exposure #{cfg.exp} ...")
    else:
        Quit(cfg, -1, f"Invalid exposure. Requesting exp #{cfg.exp} but {cfg.datevshot} has only {cfg.numexp}")

rc = initial_setup(cfg)

if rc < 0:
    Quit(cfg,rc,"Could not complete initial setup.")



#########
# after the initial setup, move stdout and stderr to a log file
#########

print(f"Starting reduction {cfg.datevshot} ... ")

#cfg.orig_stdout = sys.stdout
#cfg.orig_stderr = sys.stderr
cfg.file_stdout = open(f"{cfg.datevshot}.log","a")
print(f"Logging redirected to: {cfg.cwd}/{cfg.file_stdout.name}")
#sys.stderr = cfg.file_stdout
#sys.stdout = cfg.file_stdout


###########
# step1
###########
if s01_run1s:
    run_run1s(cfg)

    # todo: run any checks
    check_run1s(cfg)

    #todo: this would be manual here, I think, but CAN copy /red1/xxx to /scratch/local/projects
    #  all the various CoFe*.fits and multi*.fits ... these are also in the
    #  local d<shot><exp> folder in the two tar files (_co.tar and _mu.tar for the CoFe*.fits and multi*.fits respectively)
else:
    print("Skipping run1s")

###########
# step2
# VDRP
###########

if s02_vdrp:
    run_vdrp(cfg)

    check_vdrp(cfg)

    #todo: optional manual step here (need to be hetdex user), copy the *.dithall to
    #  /scratch/projects/hetdex/detect/dithall   (and /coral-repl/...)

else:
    print("Skipping vdrp")



###########
# step3
# detect (Flux Calibration)
# rallcal stuff
###########

if s03_fluxcal:

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
        print(f"flux calibration: {star_cat}")
        run_fluxcalibration(cfg,star_cat)

        if check_fluxcalibration(cfg) < 0:
            print(f"flux calibration: {star_cat} failed. Trying next ... ")
        else:
            break #this one was good

    #todo: optional:  copy detect/tp/yyyymmddvssssedtp_f.dat  to /scratch/projects/hetdex/detectp and /corral-repl/xxx

    #todo: optional: update  /scratch/projects/hetdex/detect/fwhm.all and norm.all
    #                see update_fwhm_norm script
    #

else:
    print("Skipping flux calibration")


###########
# step4
# Sky Subtractions, flux calibrations
# getcen and more stuff
###########

if s04_sky_subtraction:

    print("Getting IFU centers ...")
    if run_make_ifucen(cfg) != 0:
        Quit(cfg,-1,"FATAL. Failed to get IFU centers.")

    if s04b_rfft:
        print("Running rfft (this may take a while) ...")
        if run_rfft(cfg) != 0:
            Quit(cfg, -1, "FATAL. rfft fail. One or more expected outputs failed.")
    else:
        print("Skipping rfft")

    if s04c_rcal_all:
        print("Running rcal_all ...")
        rc = run_rcal(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rcal_all fail.")
        elif rc > 0:
            print("rcal_all: Limited success. Non-fatal. Will continue")
        #else keep going

    else:
        print("Skipping rcal_all")

    if s04d_shot_h5:
        rc = build_shot_h5(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. Could not build shot h5 file. Cannot continue with catalog creation.")

    # check stats
    if s04e_amp_stats:
        rc = amp_stats(cfg)
        if rc < 0:
            print("Non-fatal. Could not compute amp stats from shot h5 file. Will continue anyway with catalog creation.")
else:
    print("Skipping sky subtraction")


###########
# step5
# line detections
# continuum detections
###########

if s05_detection:

    if s05b_rdet_rf1:
        print("Running rdet_rf1 (line detection) ...")
        rc = rdet_rf1(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rdet_rf1 fail.")
        elif rc > 0:
            print("rdet_rf1: Limited success. Non-fatal. Will continue")
    else:
        print("skipping rdet_rf1 (line detection)")

    if s05c_rgetmax:
        print("Running rgetmax (continuum detection) ...")
        rc = rgetmax(cfg)
        if rc < 0:
            Quit(cfg, -1, "FATAL. rgetmax fail.")
        elif rc > 0:
            print("rgetmax: Limited success. Non-fatal. Will continue")
    else:
        print("skipping rgetmax (continuum detection)")


    if s05d_detection_tables:
        rc = build_detection_tables(cfg)

    if s05e_detection_hdf5:
        rc = build_detection_hdf5(cfg)


else:
    print("Skipping detections")


##################################################
# step6
# combine detections
# run diagnose and elixer
# make source catalogs
#################################################

if s06_catalogs:

    print("Catalog creation ... ")

    try:
        if s06b_fof:
            #line_tab = Table.read(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.fits"),format="fits")
            line_h5 = tables.open_file(os.path.join(cfg.cwd,f"{cfg.datevshot}_line.h5"))
            line_tab = Table(line_h5.root.Detections.read())
            line_h5.close()
            if line_tab is not None:
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
                        print(f"Updated raw lines table: {os.getcwd()}/{tname}")

                        esel = np.array(line_tab['sel_det']==True)
                        esel = esel & np.array(line_tab['continuum'] >= -3)
                        esel = esel & np.array(line_tab['sn'] >= 4.8)
                        esel = esel & np.array(line_tab['chi2'] <= 2.5)
                        #this is a bit more liberal than standard HETDEX_API
                        esel = esel & np.array(line_tab['linewidth'] >= 1.2) & np.array(line_tab['linewidth'] <= 16)

                        np.savetxt('elixer_line.dets',line_tab['detectid'][esel],fmt="%d")

                    else:
                        print("Error! (1) Could not combine lines detections by FoF.")
                else:
                    print("Error! (2) Could not combine lines detections by FoF.")
            else:
                print("Error! Could not combine lines detections. Lines table not found.")
    except:
        print("Error! Could not combine lines detections by FoF.")
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
                                print(f"Error! Unexpected matches ({np.count_nonzero(sel)}) for {cont_tab['detectid'][i]}")
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
                    print(f"Updated raw continuum table: {os.getcwd()}/{tname}")
                    #lastly sub select based on some minimum contflux ??

                    esel = np.array(cont_tab['sel_det'] == True)
                    np.savetxt('elixer_cont.dets', cont_tab['detectid'][esel],fmt="%d")


                else:
                    print("Error! (1) Could not combine continuum detections by FoF.")
            else:
                print("Error! (2) Could not combine continuum detections by FoF.")
        else:
            print("Error! Could not combine continuum detections. Continuum table not found.")

    except:
        print("Error! Could not combine continuum detections by FoF.")
        print(traceback.format_exc())


    #we now have clustered lines and continuum
    #the additional clustering performed in hetdex_api is overkill or unnecessary here, I think, as this is just
    # a single shot where hetdex_api is clustering over all the shots
    #It might even not really be necessary for the 3D and 2D clustering calls ... I think the 3D is sufficient, but
    #  for now, leave them both in and combine

    #all we REALLY care about at this point is the detectids with matching source_ids
    #we want to roll those in with the original tables and then, for each source_id just
    #  use the one with the highest SNR or lineflux? for lines and the highest continuum for cont sources



##########
# DONE
##########

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
