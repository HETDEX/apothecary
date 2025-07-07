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

import traceback

########################################################################
# CONFIGURATION
########################################################################
EchoCmds = False #if True echo system commands to the log
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot

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


#execute steps
s01_run1s = False
s02_vdrp = False
s03_fluxcal = True


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
        pass

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
                print(f"cp2dithall ... fail. File does not exist {dithall_use}")
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
        print("cp2dithall ... fail")
    else:
        print("cp2dithall ... OK")
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
cfg.file_stdout = open(f"{cfg.datevshot}.log","w")
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
