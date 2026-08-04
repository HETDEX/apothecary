#!/usr/bin/env python3
"""
helper for ds9 to display a specific frame
"""

import sys
#import glob
import os


args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

if "multiframe" in args:
    print("-multiframe switch incompatible. Do Not Use.")
    exit(0)

if "-frame" in args:
    i = args.index("-frame")
    try:
        frame = args[i+1]
    except:
        print(f"Invalid --frame specified")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-frame")
else:
    frame = None

#print(f"remaining args: {' '.join(args)}")

fns = []
del_idx = []
for i,arg in enumerate(args):
    if ".fits" in arg: #assume this is a file
        fns.append(arg)
        del_idx.append(i)

#now clean up , deleta from the end
for i in del_idx[::-1]:
    del args[i]

#print(f"searching: {args[0]}")
#fns = glob.glob(args[0])
#if len(fns) == 0:
#    print("None found")
#    exit(0)

# print(f"adding frame [{args[1]}]")
# frame = args[1]


#print(f"remaining args: {' '.join(args)}")

cmdstr = "ds9 "
for arg in args:
    cmdstr +=f" {arg}"

if frame is not None:
    for fn in fns:
        cmdstr += f" {fn}[{frame}]"
else:
    for fn in fns:
        cmdstr += f" {fn}"

#print(cmdstr)
#print("...")
os.system(cmdstr)