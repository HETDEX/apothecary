import sys
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import tables as tb
from astropy.table import Table

from hetdex_api.shot import open_shot_file, get_fibers_table

HDRVersion = "hdr5"

if __name__ == '__main__':
    
    shotid = int(sys.argv[1])
    
    try:
#        with h5py.File(f"/scratch/05865/maja_n/ffskysub/improved_spectra/mask/{shotid}.h5", "r") as ff:
        with h5py.File(f"./output/ffskysub/improved_spectra/mask/{shotid}.h5", "r") as ff:
            print(ff['calfib_ffsky_rescor'].shape)
            improved_spec = ff['calfib_ffsky_rescor'][:]
    except:
        improved_spec = None

    try:
#        with h5py.File(f"/scratch/05865/maja_n/ffskysub/improved_spectra/nomask/{shotid}.h5", "r") as ff:
        with h5py.File(f"./output/ffskysub/improved_spectra/nomask/{shotid}.h5", "r") as ff:
            print(ff['calfib_ffsky_rescor'].shape)
            improved_spec_nomask = ff['calfib_ffsky_rescor'][:]
    except:
        improved_spec_nomask = None

    try:
#        fileh = tb.open_file(f"/scratch/05865/maja_n/ffskysub/improved_spectra/pytables/nomask/{shotid}.h5", "r")
        fileh = tb.open_file(f"./output/ffskysub/improved_spectra/pytables/nomask/{shotid}.h5", "r")
        tab = Table(fileh.root.calfib_ffsky_rescor.read())
        improved_spec_nomask = tab['calfib_ffsky_rescor']
        fileh.close()
    except:
        improved_spec_nomask = None
    
    if improved_spec_nomask is None:
        N = 3
    elif improved_spec is None:
        N = 3
        improved_spec = improved_spec_nomask
        improved_spec_nomask = None
    else:
        N = 4

    N += 1 #add one more for the difference (always)

    fibers = get_fibers_table(shotid,survey=HDRVersion)
    fibers['calfib_ffsky'][fibers['calfib_ffsky']==0] = np.nan
    fibers['calfib'][fibers['calfib']==0] = np.nan

    fig = plt.figure(figsize=(17, 300))

    n = 0

    if True: #3rd column is the local sky (calfib)
        n +=1
        ax0 = fig.add_subplot(1,N,n)
        ax0.imshow(fibers['calfib'], vmin=-0.4, vmax=0.4)
        ax0.set_title("calfib")

        #fig.tight_layout()
        #fig.savefig(f"./output/plots/{shotid}.png", dpi=100, bbox_inches='tight')

    n += 1
    ax1 = fig.add_subplot(1,N,n)
    ax1.imshow(fibers['calfib_ffsky'], vmin=-0.4, vmax=0.4)
    ax1.set_title("calfib_ffsky")
    indices = np.arange(0, len(fibers['calfib_ffsky']), 112)
    ax1.set_yticks(indices, [x.replace("multi_", "") for x in fibers['multiframe'][indices]], rotation=60)

    n += 1
    ax2 = fig.add_subplot(1,N,n)
    ax2.imshow(improved_spec, vmin=-0.4, vmax=0.4)
    ax2.set_title("After residual correction.")

    if improved_spec_nomask is not None:
        n += 1
        ax4 = fig.add_subplot(1,N,n)
        ax4.imshow(improved_spec_nomask, vmin=-0.4, vmax=0.4)
        ax4.set_title("After residual correction (no mask)")


    # if True: #3rd column is the local sky (calfib)
    #     ax3 = fig.add_subplot(1,N,N)
    #     ax3.imshow(fibers['calfib'], vmin=-0.4, vmax=0.4)
    #     ax3.set_title("calfib")
    #
    #     fig.tight_layout()
    #     fig.savefig(f"./output/plots/{shotid}.png", dpi=100, bbox_inches='tight')
    if True: #3rd column is the difference between the original ffsky and the rescor ffsky
        n += 1
        ax3 = fig.add_subplot(1,N,n)

        ffsky_diff = np.nan_to_num(fibers['calfib_ffsky']) - np.nan_to_num(improved_spec)
        mx_diff = np.max(abs(ffsky_diff))
        num_diff = np.count_nonzero(ffsky_diff)
        print(f"Maximum |difference|: {mx_diff}")
        print(f"Number of different elements: {num_diff}")

        ax3.imshow(ffsky_diff, vmin=-0.4, vmax=0.4)
        ax3.set_title(f"delta, mx ({mx_diff}), #  {num_diff}")

        fig.tight_layout()

        if not os.path.exists("./output/plots"):
            #print(f"mkdir {save_filename_dir}")
            os.makedirs("./output/plots")

        fig.savefig(f"./output/plots/{shotid}d.png", dpi = 100, bbox_inches='tight')
