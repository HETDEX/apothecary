import glob
import numpy as np
from astropy.table import Table,vstack

fnpattern = "x??_source_catalog_5.0.2.fits"
merge_name = "source_catalog_5.0.3.fits"
fns = glob.glob(fnpattern)

print(f"Merging {len(fns)}...")
T = vstack([Table.read(fn) for fn in fns])
T['source_id'] = T['source_id'] + np.int64(1e10)
print(f"Writing {merge_name} ...")
T.write(merge_name,format="fits")
