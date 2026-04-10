import glob
import numpy as np
from astropy.table import Table,vstack

fnpattern = "x??_source_catalog_5.0.2.fits"
merge_basename = "source_catalog_5.0.3"
fns = glob.glob(fnpattern)

print(f"Merging {len(fns)}...")
T = vstack([Table.read(fn) for fn in fns])
T['source_id'] = T['source_id'] + np.int64(1e10)
T.add_index('detectid')
T.add_index('source_id')

print(f"Writing {merge_basename}.fits ...")
T.write(f"{merge_basename}.fits",format="fits")

#?? do we even care about an h5 version?
print(f"Writing {merge_basename}.h5 ...")
T.write(f"{merge_basename}.h5",path="SourceCatalog",format="hdf5",serialize_meta=True)


#!!!! NOTICE: this will fail if it is large and too little memory
#!!!! the ascii version requires 2x (or more) the memory of the fits
# print(f"Writing {merge_basename}.tab ...")
# T.write(f"{merge_basename}.tab",format="ascii")