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


########################################################################
# CONFIGURATION
########################################################################
FutureShotDateLimit = 20490101000  # do not allow shots after this dave+shot

ScriptRepo = "/work/03261/polonius/hetdex/single_shot"

WorkDirRoot = "./"



########################################################################
# !!! DO NOT MODIFY BELOW
#     unless you REALLY known
#     what you are doing !!!
########################################################################

@dataclass
class Config:
    clean: bool
    shotid: int
    datevshot: str
    cwd_orig: str
    cwd: str


########################################################################
# Basic user input
########################################################################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
args.pop(0) #remove THIS file

cfg = Config(False,0,"",os.getcwd(),os.getcwd())

if "-clean" in args:
    cfg.clean = True
    #idx = args.index("-clean")
    args.remove("-clean")

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

def initial_setup(cfg):
    """
    copy from script repo(s)

    change the cwd to the new workdir for this shot

    :return:
    """

    workdir = os.path.join(WorkDirRoot,f"sci{cfg.datevshot}")

    if os.path.exists(workdir):
        print("Shot directory already exists here!")
        return -1

    os.makedirs(workdir)
    os.chdir(workdir)
    cfg.cwd = workdir #now under the sci<shot> directory

    shutil.copy2(os.path.join(ScriptRepo, "science_reductions", "rsetups"),".")
    shutil.copy2(os.path.join(ScriptRepo, "science_reductions", "rfixspec"),".")


    print(f"Working under: {os.getcwd()}")
    return 0


########################################################################
# Main (section)
#   notice: no actual main function
########################################################################

rc = initial_setup(cfg)

if rc < 0:
    exit(rc)

print(f"Working under: {os.getcwd()}")



