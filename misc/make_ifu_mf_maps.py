#iterate over all shots (but need only one per date), 
#record date,  unique combos of specid, ifuslot, ifuid


#todo: add option to start from some date and append to existing ifu_mf_map.dat or _full.dat rather than overwrite

import tables
import os.path as op
import glob
import numpy as np
from astropy.table import Table
from tqdm import tqdm

#full meaning one per DAY, rather than one in month
#the assumption being that an IFU won't change over the course of 1 day, but might during the month
ifu_map_full = Table(dtype=[('yyyymmdd',np.int32),('yyyymm',np.int32),('specid', 'U3'), ('ifuslot', 'U3'),('ifuid', 'U3')])
#ifu_map = Table(dtype=[('yyyymm',np.int32),('specid', 'U3'), ('ifuslot', 'U3'),('ifuid', 'U3')])

shotdir = "/corral-repl/utexas/Hobby-Eberly-Telesco/hdr5/reduction/data"
fns = np.array(sorted(glob.glob(f"{shotdir}/2024*.h5")))

#reduce to just one per date
fns_date = [x.split("/")[-1][0:8] for x in fns]
fns_date_unique, fns_idx = np.unique(fns_date, return_index=True)
fns_unique = [op.basename(x) for x in fns[fns_idx]]
#fns = fns[fns_idx]

#one per day, but many per month
for fn in tqdm(fns[fns_idx]):
    h5 = tables.open_file(fn)
    mf = h5.root.Data.FiberIndex.read(field="multiframe")
    mf = np.unique([x.decode()[6:17] for x in mf])
    h5.close()
    
    yyyymmdd = int(op.basename(fn)[0:8])
    yyyymm = int(op.basename(fn)[0:6])
    
    for e in mf:
        toks = e.split("_")
        ifu_map_full.add_row([yyyymmdd,yyyymm,str(toks[0]).zfill(3),str(toks[1]).zfill(3),str(toks[2]).zfill(3)])
    
    # if op.basename(fn) in fns_unique:
    #     for e in mf:
    #         toks = e.split("_")
    #         ifu_map.add_row([yyyymm,str(toks[0]).zfill(3),str(toks[1]).zfill(3),str(toks[2]).zfill(3)])
    #
    #     #also write per month as we go
    #     ifu_map.write("ifu_mf_map.dat", format="ascii", overwrite=True)
    #     ifu_map_full.write("ifu_mf_map_full.dat", format="ascii", overwrite=True)

#ifu_map.add_index("yyyymm")
ifu_map_full.add_index("yyyymm")

#ifu_map.write("ifu_mf_map.dat",format="ascii",overwrite=True)
ifu_map_full.write("ifu_mf_map_full.dat",format="ascii",overwrite=True)

