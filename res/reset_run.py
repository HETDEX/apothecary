import glob
import os

try:
  run_files = glob.glob("multi*.run")
  for fn in run_files:
    open_fn = fn.split(".")[0] + ".open"
    os.rename(fn, open_fn)
except Exception as E:
  print("failed reset_run ", E)
