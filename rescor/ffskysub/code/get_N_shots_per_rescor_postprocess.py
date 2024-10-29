import sys
import os
import glob
import h5py
import numpy as np
import multiprocessing as mp

month = sys.argv[1]

#combined_pattern = f"/scratch/05865/maja_n/ffskysub/rescor/nomask/{month}/combined/combined_multi_???_???_???_??.h5"
#individual_pattern = "/scratch/05865/maja_n/ffskysub/rescor/nomask/{}/individual/*_{}_?.h5"

combined_pattern = f"/scratch/03261/polonius/rescor/code/rescor/ffskysub/rescor/nomask/{month}/combined/combined_multi_???_???_???_??.h5"
individual_pattern = "/scratch/03261/polonius/rescor/code/rescor/ffskysub/rescor/nomask/{}/individual/*_{}_?.h5"

read_filenames = np.sort(glob.glob(combined_pattern))
multiframes = [x.split("/")[-1][9:-3] for x in read_filenames]

print(f"{combined_pattern = }")

def get_N_files(i):
    ind_filenames = glob.glob(individual_pattern.format(month, multiframes[i]))
    N_files_here = len(ind_filenames)

    if i % 10 == 0:
        print(f"Finished {i}/{N}")
    return N_files_here

N = len(read_filenames)
with mp.Pool(8) as pool:
    N_files = pool.map(get_N_files, range(N))

outfilename = os.path.join("/".join(combined_pattern.split("/")[:-1]), "N_files_postprocess.txt")
with open(outfilename, "w") as ff:
    for i in range(len(multiframes)):
        ff.write(f"{multiframes[i]} {N_files[i]}\n")
print(f"Wrote to {outfilename}")
