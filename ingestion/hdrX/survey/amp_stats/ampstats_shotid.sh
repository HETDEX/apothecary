#!/usr/bin/env python

from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import sys
import time
import os.path as op


start = time.time()

cl_args = list(map(str.lower,sys.argv))

shotid = int(cl_args[1])

if op.exists(f"{shotid}_stats.pickle"):
    print(f"[{shotid}] already done. Skipping.")
    exit(0)

print(f"[{shotid}] laoding stats, start: {time.strftime('%H:%M:%S')}")
_ = AS.make_stats_for_shot(shotid=shotid,preload=True)
stop = time.time()
print(f"[{shotid}] elapsed time: {stop-start}")
exit(0)

try:
  sd = AS.load_shot_stats_pickle(shotid) #if this actually loads, the pickle exists and this one was already done
  if sd is None or len(sd) == 0:
    print(f"[{shotid}] starting: {time.strftime('%H:%M:%S')}")
    print(f"[{shotid}] loading stats (1) ....")
    _ = AS.make_stats_for_shot(shotid=shotid,preload=True)
  else:
    print(f"[{shotid}] already done. Skipping.")
except:
  print(f"[{shotid}] loading stats (2) ....")
  _ = AS.make_stats_for_shot(shotid=shotid,preload=True)

stop = time.time()
print(f"[{shotid}] elapsed time: {stop-start}")
