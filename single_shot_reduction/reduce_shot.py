"""
Perform science reduction (astrometry, flux calibration, line and continuum detection, elixer )
  on a single shot (field). Can be HETDEX dither-style or not.

  Copies code source and works in the current directory

  input: [required] shotid (or datevshot)
         [optional] -clean  (clean up workfiles, leaving only the output and logs)



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


DEBUG_NO_COPY = False

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



########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:
    clean: bool
    overwrite: bool
    shotid: int
    datevshot: str
    cwd_orig: str
    cwd: str
    scriptdir: str


########################################################################
# Basic user input
########################################################################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

cfg = Config(False,False,0,"",os.getcwd(),os.getcwd(),"")

if "-clean" in args:
    cfg.clean = True
    #idx = args.index("-clean")
    args.remove("-clean")

if "-overwrite" in args:
    cfg.overwrite = True
    args.remove("-overwrite")

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
    if msg is not None:
        print(f"({rc})",msg)
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
            if DEBUG_NO_COPY:
                print("DEBUG ... skip overwrite directory")
            else:
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
        cfg.scriptdir = ScriptRepo

    os.chdir(workdir)
    cfg.cwd = workdir #now under the sci<shot> directory

    if DEBUG_NO_COPY:
        print("DEBUG ... skip source code copy")
        return 0

    print("Copying source code ...")

    #shutil.copy2(os.path.join(cfg.scriptdir, "science_reductions", "rsetups"),".") #no, this function is its equivalent

    shutil.copy2(os.path.join(cfg.scriptdir, "rfixspec"), ".")
    shutil.copytree(os.path.join(cfg.scriptdir, "sciscripts"), ".", dirs_exist_ok=True)
    shutil.copytree(os.path.join(cfg.scriptdir,"vdrp"), "vdrp", dirs_exist_ok=True)



    return 0


def run1s():
    """

    science_reduction/sciscripts/run1  batch file


    :return:
    """




########################################################################
# Main (section)
#   notice: no actual main function
########################################################################

rc = initial_setup(cfg)

if rc < 0:
    Quit(cfg,rc,"Could not complete initial setup.")





