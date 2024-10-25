#!/usr/bin/env python

from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
import sys


#usually a YYYYMM
cl_args = list(map(str.lower,sys.argv))

try:
  fn = op.basename(cl_args[1])
  shotid = int(fn[0:8] + fn[9:12])
  sd = AS.load_shot_stats_pickle(shotid) #if can load the pickle, it already has been done, so skip
  if sd is None or len(sd) ==0:
    print(f"{shotid} pickle not found, will make stats")
    try:
      _ = AS.make_stats_for_shot(fqfn=cl_args[1],preload=True)
    except:
      print(f"Exception in make_stats_for_shot")
  else:
    print(f"{shotid} already done. Skipping.")
except:
  print(f"{shotid} pickle not found, will make stats")
  try:
    _ = AS.make_stats_for_shot(fqfn=cl_args[1],preload=True)
  except:
    print(f"Exception in make_stats_for_shot")
