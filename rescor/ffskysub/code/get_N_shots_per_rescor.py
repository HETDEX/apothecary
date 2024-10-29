import sys
import os
import glob
import h5py
import numpy as np

month = sys.argv[1]

read_pattern = f"/scratch/05865/maja_n/ffskysub/rescor/nomask/{month}/combined/combined_multi_???_???_???_??.h5"
read_filenames = np.sort(glob.glob(read_pattern))
multiframes = [x.split("/")[-1][9:-3] for x in read_filenames]

N_files = []
for read_filename in read_filenames:
    with h5py.File(read_filenames, "r") as ff:
        N_files.append(ff['rescor'].attrs['N_files'])

outfilename = os.path.join("/".join(read_pattern.split("/")[:-1]), "N_files.txt")
with open(outfilename, "w") as ff:
    for i in range(len(multiframes)):
        ff.write(f"{multiframes[i]} {N_files[i]}\n")
print(f"Wrote to {outfilename}")