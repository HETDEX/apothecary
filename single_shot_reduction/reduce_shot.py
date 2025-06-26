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
import shutil
from dataclasses import dataclass

import traceback

########################################################################
# CONFIGURATION
########################################################################
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot

ScriptRepo = "/work/03261/polonius/hetdex/single_shot"
LocalScriptRepo = "./local_script_repo" #useful if running multiple single shots ... can copy remotely once
                                        #then copy locally from here for each shot
                                        #set to None if you do NOT want to use a local script dir cache
                                        #  and force a copy from the main repo each time

WorkDirRoot = "./"

HETRaw = "/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw/"
karlgettar = "/work/00115/gebhardt/maverick/gettar/"
karlhome = "/home1/00115/gebhardt"

########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:

    clean: bool = False
    overwrite: bool = False
    shotid: int = 0
    datevshot: str = ""
    exp: int = 0  #specific exposure number to reduce
    numexp: int = 0 #number of exposures in the shot
    cwd_orig: str = os.getcwd()
    cwd: str = os.getcwd()
    scriptdir: str = ""
    gettar_fn: str =  ""  # the runs* or runt* file from karlgettar folder with the date, shot, exp data

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

    if cfg.orig_stdout:

        sys.stdout = cfg.orig_stdout

    if cfg.orig_stderr:
        sys.stderr = cfg.orig_stderr

    if cfg.file_stdout:
        cfg.file_stdout.close()

        #repeat the final message to the console
        if msg is not None:
            print(f"({rc})", msg)
        else:
            print(f"({rc})")

    exit(rc)


def initial_setup(cfg):
    """
    copy from script repo(s)

    change the cwd to the new workdir for this shot

    this is largely equivalent to what rsetups would do, but here for a single shot and not a whole month

    :return:
    """

    workdir = os.path.join(WorkDirRoot,f"sci{cfg.datevshot}")

    if os.path.exists(workdir):
        if cfg.overwrite:
            print(f"Overwriting directory {workdir} ... ")
            shutil.rmtree(workdir)
        else:
            print(f"Shot directory already exists here! {workdir}")
            return -1

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

    os.chdir(workdir)
    cfg.cwd = os.getcwd() #now under the sci<shot> directory


    print("Copying source code ...")

    #shutil.copy2(os.path.join(cfg.scriptdir, "science_reductions", "rsetups"),".") #no, this function is its equivalent

    shutil.copy2(os.path.join(cfg.scriptdir, "rfixspec"), ".")
    shutil.copytree(os.path.join(cfg.scriptdir, "sciscripts"), ".", dirs_exist_ok=True)
    shutil.copytree(os.path.join(cfg.scriptdir,"vdrp"), "vdrp", dirs_exist_ok=True)

    #update the "home" path tilde
    os.system(f"sed -i s#~gebhardt#{karlhome}# rbfits")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rbfits_fix")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rback_field")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rback_fix")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rbacks")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rbfits_s")  # use '#' as sed separator rather than "/"
    os.system(f"sed -i s#~gebhardt#{karlhome}# rimarb")  # use '#' as sed separator rather than "/"

    #os.system(f"sed -i s/ChangeMe/${mth}/ run1s")

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


def run1s(cfg):
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
        #run1s 20240730 009 exp01 202407

        #sed -i s/\.\.\/runsChangeMe/${mth}/ run1s
        #sed -i s/ChangeMe/${mth}/ run2s
        #sed -i s/ChangeMe/${mth}/ rtaremc


        os.system(f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# run1s") #use '#' as sed separator rather than "/"
        os.system(f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# run2s")  # use '#' as sed separator rather than "/"
        #os.system(f"sed -i s#../runsChangeMe#{cfg.gettar_fn}# rtaremc")  # use '#' as sed separator rather than "/"

        #cmd = "sed -i s#\${scriptdir}"+f"#{cfg.scriptdir}/sciscripts/# run1s"
        #scripts have already been copied to shot workding dir
        cmd = "sed -i s#\${scriptdir}" + f"#{cfg.cwd}/# run1s"
        os.system(cmd)  # use '#' as sed separator rather than "/"

        #actually run it here
        os.system(f"run1s {cfg.datevshot[0:8]} {cfg.datevshot[-3:]} exp{str(exp).zfill(2)} {cfg.datevshot[0:6]}")





########################################################################
# Main (section)
#   notice: no actual main function
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

cfg.orig_stdout = sys.stdout
cfg.orig_stderr = sys.stderr
cfg.file_stdout = open(f"{cfg.datevshot}.log","w")
sys.stderr = cfg.file_stdout
sys.stdout = cfg.file_stdout



print("this is just a log test ....")

###########
# step1
###########
#run1s(cfg)






##########
# DONE
##########

Quit(cfg, 0, f"Complete: {cfg.datevshot}")