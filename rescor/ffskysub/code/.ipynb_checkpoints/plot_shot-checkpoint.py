import sys
import h5py
import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    shotid = int(sys.argv[1])
    
    with h5py.File(f"/scratch/05865/maja_n/ffskysub/improved_spectra/{shotid}.h5", "r") as ff:
        print(ff['calfib_ffsky_rescor'].shape)
        improved_spec = ff['calfib_ffsky_rescor'][:]
        
    fibers = get_fibers_table(20230826016)
    fibers['calfib_ffsky'][fibers['calfib_ffsky']==0] = np.nan
    fibers['calfib'][fibers['calfib']==0] = np.nan

    fig = plt.figure(figsize=(20, 300))

    ax1 = fig.add_subplot(131)
    ax1.imshow(fibers['calfib_ffsky'], vmin=-0.4, vmax=0.4)
    ax1.set_title("calfib_ffsky")
    indices = np.arange(0, len(fibers['calfib_ffsky']), 112)
    ax1.set_yticks(indices, [x.replace("multi_", "") for x in fibers['multiframe'][indices]], rotation=60)

    ax2 = fig.add_subplot(132)
    ax2.imshow(improved_spec, vmin=-0.4, vmax=0.4)
    ax2.set_title("After residual correction.")

    ax3 = fig.add_subplot(133)
    ax3.imshow(fibers['calfib'], vmin=-0.4, vmax=0.4)
    ax3.set_title("calfib")

    fig.tight_layout()
    fig.save_fig(f"../plots/{shotid}.png", dpi = 400, bbox_inches='tight')