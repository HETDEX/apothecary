"""
usage: python sim_mc2txt.py <shotid>
"""


import tarfile as tar
from tqdm  import tqdm
import numpy as np
import copy
import os.path as op
from astropy.table import Table
import sys

if op.exists("/home/jovyan/Hobby-Eberly-Telesco/"):
    tarbasedir = "/home/jovyan/Hobby-Eberly-Telesco/simulations/xflim_archive/"
else:
    tarbasedir = "/corral/utexas/Hobby-Eberly-Telesco/simulations/xflim_archive/"

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file

try:
    shotid = int(args[0]) #202408
except:
    shotid = int(args[0].replace('v','').replace('s','').replace('d',''))


DTYPE = np.dtype([
    ("line_num", "<i4"), ("wave", "<i4"), ("xflim_bg", "<f4"),
    ("Fout", "<f4"), ("LWout", "<f4"), ("SN", "<f4"), ("xflim", "<f4"),
])

batch_ids = np.arange(0,68)

# find the right tar file
shotstr = f"{str(shotid)[0:8]}v{str(shotid)[-3:]}".encode()
matched_batch_id = -1
batch_tarfile = None

print(f"Searching for the batch with shotid: {shotid} ...")
with tar.open(op.join(tarbasedir, "xflim_pkg.tar"), "r") as tarfh:
    # all_names = tarfh.getnames()
    for bn in batch_ids:
        batch = f"xflim_pkg/batches/batch_{str(bn).zfill(2)}.list"
        with tarfh.extractfile(batch) as f:
            if shotstr in f.read():
                matched_batch_id = bn
        if matched_batch_id > -1:
            break

if matched_batch_id > -1:
    print("Found")
    batch_tarfile = op.join(tarbasedir, f"xflim_batch_{str(matched_batch_id).zfill(2)}.tar")

    print(f"Extracting and converting mc files from: {batch_tarfile} ...")
    with tar.open(batch_tarfile, "r") as tarfh:
        shot_path_in_tar = f"sim{str(shotid)[0:8]}v{str(shotid)[-3:]}s666anchor/output"
        mc_files = [name for name in tarfh.getnames() if name.startswith(shot_path_in_tar)]
        for mc in tqdm(mc_files):
            if mc[-3:] == ".mc":
                with tarfh.extractfile(mc) as f:
                    x = copy.copy(np.frombuffer(f.read(), dtype=DTYPE))
                    # convert wave_idx to wave
                    x['wave'] = 3470 + (x['wave'] - 1) * 2
                    # writeout
                    outfn = op.basename(mc).replace(".mc", ".txt")
                    Table(x).write(outfn, format="ascii", overwrite=True)
else:
    print("Not found")