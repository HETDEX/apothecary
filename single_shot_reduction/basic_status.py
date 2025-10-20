"""
 walk the sci* directories, open the shot h5 file and spit out basic status

 shotid Ra dec total_time #exp seeing transmission

"""


import glob
import tables
import traceback
import numpy as np
import os


h5_fns = sorted(glob.glob("sci*/20*v???.h5"))
print(f"#shotid\t\tra\t\tdec\t\texp\ttime\tfwhm\ttput\tline_ct\tcont_ct\tdither")
for fn in h5_fns:
    cont_ct = -1
    line_ct = -1

    try:
        sci = fn.split("/")[0]

        det_fn = os.path.join(sci,"elixer/line.dets")
        if os.path.exists(det_fn):
            with open(det_fn,'r') as f:
                lines = f.read()
                line_ct = lines.count("\n")

        det_fn = os.path.join(sci,"elixer/cont.dets")
        if os.path.exists(det_fn):
            with open(det_fn,'r') as f:
                lines = f.read()
                cont_ct = lines.count("\n")

    except:
        print(f"[{fn}] Exception (trying wc):", traceback.format_exc())

    try:
        h5 = tables.open_file(fn,mode="r")
        row = h5.root.Shot.read()[0]  #there should only ever be one row

        times = row['exptime']
        num_exp = np.count_nonzero(times)
        total_time = np.sum(times)

        xs = [f"{x:0.2f}" for x in row['xditherpos']]
        ys = [f"{y:0.2f}" for y in row['yditherpos']]
        dither = [(x,y) for x,y in zip(xs,ys)]

        print(f"{row['shotid']}\t{row['ra']:0.6f}\t{row['dec']:0.6f}\t{num_exp}\t{total_time:0.1f}\t"
              f"{row['fwhm_virus']:0.2f}\t{row['response_4540']:0.2f}\t{line_ct}\t{cont_ct}\t{dither}")

        h5.close()
    except:
        print(f"[{fn}] Exception:", traceback.format_exc())


