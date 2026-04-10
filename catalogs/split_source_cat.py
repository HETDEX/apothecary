
from astropy.table import Table
from hetdex_api.extinction import *  #includes deredden_spectra also gets config?
from tqdm.notebook import tqdm


#basedir = "/home/jovyan/Hobby-Eberly-Telesco/hdr5/catalogs"
basedir = "/scratch/projects/hetdex/hdr5/catalogs"
srccat_name = "source_catalog_5.0.2.fits"
srcat_format = "fits"
num_chunks = 50
survey="hdr5"

SrcTab = Table.read(f"{basedir}/{srccat_name}")#load_source_table(flag_best=False) #Note, will apply flag_best later
SrcTab.add_index('detectid')
SrcTab.add_index('source_id')
print(f"Num detections: {len(SrcTab)}")


#get the unique source_ids
u_source_ids = np.unique(list(SrcTab['source_id']))
print(f"Num sources: {len(u_source_ids)}")

print(f"Breaking into {num_chunks} chunks ...")
src_list = np.array_split(u_source_ids,num_chunks)
t = 0
for i,s in enumerate(src_list):
    print(i, len(s))
    t = t + len(s)
#print(t)

#print(f"Writing out ...")
#make sub-source_catalogs
for i in tqdm(range(len(src_list))):
    sel = np.in1d(SrcTab['source_id'],src_list[i])

    SubSrcTab = SrcTab[sel]

    #sanity check
    if sorted(np.unique(SubSrcTab['source_id'])) != sorted(src_list[i]): #this is a problem
        print(f"Fail. Sublists do not match. index [{i}]")
        break

    #write out SubSrcTab with incrementing name
    SubSrcTab.write(f"x{str(i).zfill(2)}_{srccat_name}",format=srcat_format)
    del SubSrcTab

print("Done")