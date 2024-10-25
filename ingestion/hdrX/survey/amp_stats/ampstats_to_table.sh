#!/usr/bin/env python

from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
import sys
import glob
from tqdm import tqdm

fns = sorted(glob.glob("*_stats.pickle"))
shot_dict_list = []

print("loading pickles ...")
for fn in tqdm(fns):

  shotid = int(op.basename(fn).split("_")[0])
  sd = AS.load_shot_stats_pickle(shotid,path="./")
  shot_dict_list.append(sd)

print("building dat file (old style) ...")
outfile = "ampstats.dat"
tb = AS.stats_save_as(shot_dict_list,outfile,oldstyle=True,header=True)
print(f"Done. {outfile}")

print("building dat file (as fits) ...")
outfile = "ampstats.fits"
tb = AS.stats_save_as(shot_dict_list,outfile,format="fits",oldstyle=False,header=False)
print(f"Done. {outfile}")



