#!/usr/bin/env python

#interactive 
#commandline has one required parameter, the filename of the table to update
#next parameters, all optional, in order: shotid, multiframe (use 0 for all), expid (use 0 for all), flag, reason
#if any is bad, all after are considered bad


from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
import sys



#cl_args = list(map(str.lower,sys.argv))
cl_args = list(sys.argv) #no, we need this to be case sensitive for the strings

#has to be there
fn = cl_args[1]

shotid = None
mf = None
exp = None
flag = None
reason = None


try:
  shotid = int(cl_args[2])
  mf = cl_args[3]
  if mf == 0 or (len(mf) == 1 and mf[0] == "0"):
    mf = int(0)

  exp = int(cl_args[4])
  if exp == 0 or exp =="0":
    exp = int(0)

  flag = int(cl_args[5])
  reason = cl_args[6]


except:
  pass
  #print("Exception case")


if shotid is None:
  shotid = input(f"shotid? (integer format): ")
  shotid = int(shotid)

if mf is None:
  mf = input(f"multiframe?: ")
  if mf == 0 or (len(mf) == 1 and mf[0] == "0"):
    mf = None
elif mf == 0:
  mf = None

if exp is None:
  exp = input(f"exposure num?: ")
  exp = int(exp)
  if exp == 0 :
    exp = None

elif exp == 0:
  exp = None

if flag is None:
  flag = input(f"manual flag value?: ")
  flag = int(flag)

if reason is None or len(reason) == 0:
  reason = input(f"reason/notes?: ")

AS.stats_update_flag_manual(fn, shotid, multiframe=mf,expnum=exp,flag_manual=flag, flag_manual_desc=reason, savefmt="fits", interactive=True)


