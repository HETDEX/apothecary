import glob
import shutil

import numpy as np
from astropy.table import Table, Column
import tables

import healpy as hp
import tarfile as tar
import os
#import io
from pathlib import Path
from astropy.table import vstack
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm
from hetdex_api.extinction import deredden_spectra
import traceback
#from elixer import utilities as Utils


##################################################
# setup
##################################################

out_tab_name = "virus_shot_summary_table.update.tab"
basepath = "/work/03946/hetdex/maverick/"
mindate = 20251001
maxdate = 20261231 #last contiguous date on corral's het_raw


#make sure have the list files:
list_src_path = "/work/03261/polonius/hub/parallel"
for listfn in ["shotlist.hetdex","shotlist.include"]:
    if not os.path.exists(f"./{listfn}"):
        try:
            shutil.copy(os.path.join(list_src_path,listfn),f"./{listfn}")
        except:
            print(f"Failed to copy {listfn} from {list_src_path}\n{traceback.format_exc()}")



####################################
# optional resume
####################################
try:
    T =  Table.read(out_tab_name,format="ascii")
    lastdate = np.max(T['date'])
except:
    T = None
    lastdate = None

if lastdate is not None:
    print(f"Resuming from {lastdate}")

##############################
# helper funcs
##############################

dexlist = None
other_list = None


def get_healpix_id(ra, dec, Nside=32768):
    """
    return single integer, healpix ID, that covers this ra, dec

    separations vary, but maximum is 6.7"
    (hp.max_pixrad(32768, degrees=True) * 3600) == 6.7288)

    :param ra: degrees (float)
    :param dec: degrees (float)
    :param Nside: default 2**15 to match HETDEX_API
    :return: single integer, healpix ID
    """

    try:
        hp_id = None
        hp_id = hp.ang2pix(Nside, ra, dec, lonlat=True)
    except:
        print(f"Exception! get_healpix_id")

    return hp_id


def get_healpix_region(ra, dec, radius=6.73, Nside=32768):
    """
    return healpix region (list of healpix IDs (integers)) that contains ra, dec out to radius
    (Nside needs to match what is being searched to generate the same IDs)

    separations vary, but maximum is 6.7"
    (hp.max_pixrad(32768, degrees=True) * 3600) == 6.7288)

    :param ra: degrees (float)
    :param dec: degrees (float)
    :param radius: arcsecs (float)
    :param Nside: default 2**15 to match HETDEX_API
    :return: list of healpix ids or None
    """

    try:
        hp_region = None
        vec = hp.ang2vec(ra, dec, lonlat=True)
        hp_region = hp.query_disc(Nside, vec, (radius / 3600.0 * np.pi / 180))
    except:
        print("Exception! get_healpix_region")

    return hp_region


def exclude_observations(T):
    """
    todo: load an excusion list
    NOTICE do not KEEP these ... this returns the matches to a list that should be excluded, removed, or marked as not to be used
    """
    # no, we will do this only as an inclusion pattern + inclusion list
    vdfi_sel = np.array(["VDFI" not in x for x in T['name']])

    sel = vdfi_sel  # or , or or, ....
    return sel


def include_observations(T):
    """
    """
    # this is not the most efficient, but it is clear to read and easy to update
    # and this is rarely run

    global dexlist, other_list

    include_sel = []

    include_sel.append(np.array(["parallel".casefold() in x.casefold() for x in T['name']]))  # case insensitive
    include_sel.append(np.array(["DEX" in x for x in T['name']]))  # case sensitive
    # note: it seems there are maybd ~ 1000 HETDEX observations NOT named this way

    if dexlist is None or other_list is None:
        try:
            dexlist = np.loadtxt("shotlist.hetdex", dtype=str)
            other_list = np.loadtxt("shotlist.include", dtype=str)
        except:
            try:
                dexlist = np.loadtxt("../shotlist.hetdex", dtype=str)
                other_list = np.loadtxt("../shotlist.include", dtype=str)
            except:
                print("Failed to load dexlist or other_list")

    include_sel.append(np.array([x in dexlist for x in T['shotid']]))

    if len(other_list) > 0:
        include_sel.append(np.array([x in other_list for x in T['shotid']]))

    sel = np.any(include_sel, axis=0)
    # print(np.shape(include_sel), np.count_nonzero(sel))
    return sel


def update_shot_status(T, basedir="./", file_list=[]):
    """
    look over all shots (ssr*h5 files) in the basedir and update the corresponding shot status
     or do the same, but use the files specified in file_list (fully qualified or relative)

    """

    if file_list is not None and len(file_list) > 0:
        fns = file_list
    else:
        fns = sorted(glob.glob(os.path.join(basedir, "ssr_*.h5")))

    if "avg_sky" not in T.columns:
        # T["avg_sky"] = -1.0
        T.add_column(Column(data=-1.0, name="avg_sky", dtype="float32"))
    if "reduction_date" not in T.columns:
        # T["reduction_date"] = 99990101
        T.add_column(Column(data=99990101, name="reduction_date", dtype="int32"))
    if "reduction_status" not in T.columns:
        #  T["reduction_status"]
        T.add_column(Column(data="unk", name="reduction_status", dtype="U4"))

    for fn in tqdm(fns):
        h5 = None
        try:
            dvs = os.path.basename(fn).replace("ssr_", "").split(".")[0]
            sel = T['shotid'] == dvs
            if np.any(sel):
                h5 = tables.open_file(fn, mode='r')
                row = h5.root.Shot.read()  # only one line

                # T['reduction_status'][sel] = "test"
                T['avg_sky'][sel] = row['avg_sky']
                T['reduction_date'][sel] = row['reduction_date']
                T['reduction_status'][sel] = row['reduction_status']
                # T['reduction_status_ext'][sel] # don't bother with the extended status
                h5.close()

            else:
                print(f"[{dvs}] unexpected. Not found in Table.")

        except:
            print(traceback.format_exc())
            if h5 is not None:
                try:
                    h5.close()
                except:
                    pass


def update_HET_shots(T, sort=False):
    """

    append new shots to the list
    if sort is True, then re-sort based on shotid

    need to look in a path ? maybe hetdex/maverick??
    need to make sure the shot to be added is not already in the table
        (and fits the inclusion rules ... so can briefly, blindly include (in memory) and then, at the end, run include_observations(T)

    ?? maybe should always re-sort (should not need to re-index as that is already set)

    """

    # todo: look at dates and add to table?

    # now make sure the inclusion rules are applied
    include_observations(T)




##############################################
##############################################
# main section
##############################################
##############################################


########################################
# update loop
########################################


fps = np.array(sorted(glob.glob(f"{basepath}/20*", recursive=False)))

if lastdate is not None or maxdate is not None or mindate is not None:
    # months = np.array([int(os.path.basename(os.path.normpath(p))) for p in fps])
    months = np.array([int(os.path.basename(os.path.normpath(f))) for f in fps])

if mindate is not None:
    sel_min = months >= mindate
else:
    sel_min = np.full(len(months), True)

if lastdate is not None:
    sel_last = months >= lastdate  # since there are multiples per date, we should pick up including the "lastdate"
else:
    sel_last = np.full(len(months), True)

if maxdate is not None:
    sel_max = months <= maxdate
else:
    sel_max = np.full(len(months), True)

fps = fps[sel_min * sel_max * sel_last]

print(f"Working on: {len(fps)} from {fps[0]} to {fps[-1]}")

# NOTE: the collection of the internal tar listing can be slow ... this can take a few minutes per date

# for p in tqdm(fps[0:1]):
for p in tqdm(fps):
    # this will end up being a series of small tables, (one date each) that we merge as we go so we can start over
    # if/when the job times out
    col_shotid = []  # datevshot ... not can fail for non-science shots with long obsids, but that is not what is included here
    col_date = []  # YYYYMMDD
    col_obsid = []  # integer
    col_name = []  # OBJECT
    col_exptime = []
    col_expnum = []
    col_ra = []
    col_dec = []
    # col_num_exp = [] #may not know this from this data
    col_dust_corr = []
    col_healpix = []

    col_pa = []  # position angle
    col_az = []  # azimuth
    col_allowed = []

    fns = np.array(sorted(glob.glob(f"{p}/virus/virus*.tar")))
    # only care about science shots
    obsid = np.array([int(os.path.basename(fn).split(".")[0][5:]) for fn in fns])
    sel = obsid < 999
    fns = fns[sel]

    # for fn in fns[0:1]:
    for fn in tqdm(fns):

        # not the same as the normal paths, so cannot use the get_shotid_from_path
        # '/work/03946/hetdex/maverick/20260201/virus/virus0000001.tar'
        toks = fn.split("/")
        ix = toks.index("virus")
        date = int(toks[ix - 1])
        obsid = int(toks[ix + 1].split(".")[0][5:])
        datevshot = str(date) + "v" + str(obsid).zfill(3)

        # it is possible we already have this one ... say, on a re-run, esp. if we want to double check the last date
        if T is not None and np.count_nonzero(T['shotid'] == datevshot) > 0:
            # print(f"[{datevshot}] already show this shotid in the table. Skipping.")
            continue

        # get the file listing from the tar
        tf = tar.open(fn)
        tarpaths = np.array(tf.getnames())

        # is this a science observation?
        obstype = os.path.basename(tarpaths[0]).split('.')[-2][-3:]
        if obstype != 'sci':
            # print(f"[{obstype}] : skipping {fn}")
            continue  # this is NOT a science observation; could be dark (drk), flat (flt), sky, cmp, zro, etc,

        # debug statement
        # print(f"[{obstype}] : summarizing {fn}")

        # how many exposures? the basepath is the same, except for the file
        exp_paths, exp_path_ix = np.unique([os.path.dirname(tp) for tp in tarpaths], return_index=True)
        num_exps = len(exp_paths)

        hdu_out_path = f"./{str(date)[:-2]}"
        if not os.path.exists(hdu_out_path):  # strip off the day
            Path(hdu_out_path).mkdir(parents=True, exist_ok=True)

        # get the basic info, the first exposure and first amp is fine
        fileh = tf.extractfile(tarpaths[0])  # open the fits

        with fits.open(fileh) as hdu:
            # get the cards
            header = hdu[0].header
            # need this to be from the tarpaths[] NOT fileh.name (which holds the outer tar file's name)
            expnum = int(tarpaths[0].split("/exp")[1][0:2])

            try:
                coord = SkyCoord(header["TRAJRA"], header["TRAJCDEC"], unit=(u.hourangle, u.deg))
            except:
                try:
                    coord = SkyCoord(header["TELRA"], header["TELDEC"], unit=(u.hourangle, u.deg))
                except:
                    continue

            if coord.ra.degree < -999. or coord.dec.degree < -999.:
                col_ra.append(-999.9)
                col_dec.append(-999.9)
                col_dust_corr.append(-999.9)
                col_healpix.append(-999)
                col_pa.append(-999.9)
                col_az.append(-999.9)
            else:
                col_ra.append(round(coord.ra.degree, 6))
                col_dec.append(round(coord.dec.degree, 6))
                try:
                    coord = SkyCoord(ra=col_ra[-1] * u.deg, dec=col_dec[-1] * u.deg)
                    dust_corr = deredden_spectra([4540.], coord)[0]
                    col_dust_corr.append(round(dust_corr, 2))
                except:
                    print("Exception in dust corr")
                    col_dust_corr.append(-999.9)
                try:
                    healpix = get_healpix_id(col_ra[-1], col_dec[-1])
                    if healpix is not None:
                        col_healpix.append(healpix)
                except:
                    col_healpix.append(-999)

                # print(col_healpix[-1])

            # change to PEXPTIME as may be more appropriate, the -8.0s is from Greg Zeimann
            col_exptime.append(round(float(header["PEXPTIME"]) - 8.0, 1))
            col_expnum.append(num_exps)

            col_shotid.append(datevshot)
            col_date.append(date)
            col_obsid.append(obsid)
            col_name.append(header['OBJECT'])

            try:
                col_pa.append(header['PARANGLE'])
            except:
                col_pa.append(-999.9)
            try:
                col_az.append(header['STRUCTAZ'])
            except:
                col_az.append(-999.9)

            col_allowed.append(0)  # assume NO, will set to YES later where appropriate
            # for future followup or re-dos
            with open(f"{hdu_out_path}/{datevshot}e{str(expnum).zfill(2)}.hdu", "w") as headerfile:
                headerfile.write(header.tostring(sep='\n'))

            # fileh.close()

            # now need any other exposures to build up exptime
            for exp_path in tarpaths[exp_path_ix]:  # this is fast
                if exp_path == tarpaths[0]:
                    continue  # already have this one

                fileh = tf.extractfile(exp_path)  # open the fits
                with fits.open(fileh) as hdu:
                    expnum = int(exp_path.split("/exp")[1][0:2])
                    header = hdu[0].header

                    try:
                        next_exp_time = round(float(header["PEXPTIME"]) - 8.0, 1)
                    except:
                        next_exp_time = round(float(header["EXPTIME"]), 1)

                    col_exptime[-1] = round(col_exptime[-1] + next_exp_time)

                    with open(f"{hdu_out_path}/{datevshot}e{str(expnum).zfill(2)}.hdu", "w") as headerfile:
                        headerfile.write(header.tostring(sep='\n'))

                fileh.close()

        tf.close()

    # after each date, merge in
    t2 = Table([col_shotid, col_name, col_date, col_obsid, col_ra, col_dec, col_pa, col_az, col_expnum, col_exptime,
                col_dust_corr, col_healpix, col_allowed],
               names=['shotid', 'name', 'date', 'obsid', 'ra', 'dec', 'PA', 'Az', 'num_exp', 'exptime',
                      'dust_corr_4540', 'healpix', 'allowed'])
    t2.sort('shotid')
    if len(t2) > 0:
        sel = include_observations(t2)  # only include from the allowed list
        t2['allowed'][sel] = 1

        if T is None:
            T = t2
        else:
            T = vstack([T, t2])

        T.write(out_tab_name, format="ascii", overwrite=True)

# final write with indexes
if T is not None:
    T.add_index(['shotid', 'ra', 'dec', 'healpix'])
    try:
        # this works for astropy version > 7.0
        T.write(out_tab_name.replace(".tab", ".fits"), format="fits", write_indices=True, overwrite=True)
    except:
        T.write(out_tab_name.replace(".tab", ".fits"), format="fits", overwrite=True)


