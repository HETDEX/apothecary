"""
re-evaluate the assigned z_hetdex and clustering of individual detections into source_ids

this is inteded to be run in SLURM against chunks of HETDEX source_catalogs

this is limited use, so only minimal error control and minimal command line support

"""
import os
import sys
import numpy as np
import tables
from astropy.table import Table
from elixer import global_config as G
from hetdex_api.extinction import *  #includes deredden_spectra also gets config?
import traceback
from tqdm import tqdm

try:
    from filelock import FileLock
except:
    print("You need to install filelock (e.g.: pip install --user filelock) ")
    exit(-1)

import warnings

# Filter by the specific message or category
print("!!!! warning. Numpy deprication warning for alltrue turned off. Need to update (py)tables")
warnings.filterwarnings("ignore", category=DeprecationWarning, message=".*alltrue.*")

#PyTables >= 3.9.0 requires Python 3.9 and blosc2  2.2.8.
#PyTables 3.8.0 is the last stable version for Python 3.8, and it is frequently forced by pip when it detects the older Python version.


##########################
# edit as needed
##########################
SHOW_TQDM = True
catchunk_format = "fits" #for use with astropy table read/write
overwrite = True # overwrite the input file after correction, if False, write out as a new file
                 # note: always writes out in the cwd, so if the input is in another directory
                 # the overwrite is mostly irrelevant and only alters the output name to avoid collision

elixer_h5_fn = "/scratch/projects/hetdex/hdr5/detect/elixer.h5"
check_counterpart_dist = True #aka update_counterpart_dist
check_z_hetdex = True
atol_strong = 0.0015  # e.g. vs emission line fit
atol_weak = 0.015  # e.g. vs an SED, photz fit
source_id_counter_fn = "source_id.counter"
source_id_counter_mutex = "source_id.lock"
global_source_id_counter = np.int64(5024e9) #max in 5.0.2 source_catalog is 5023000273239

#####################
#get ths input
#####################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

if "-catchunk" in args: #path to the shot h5 file
    i = args.index("-catchunk")
    try:
        catchunk = args[i+1]
    except:
        print(f"Invalid -catchunk specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-catchunk")
else:
    print("Fatal. Must supply --catchunk <path to the cat chunk file>")
    exit(-1)


########################
# load the table
########################

try:
    SrcTab = Table.read(catchunk, format=catchunk_format)
except:
    print(f"Fatal. Could not load {catchunk} as {catchunk_format}\n",traceback.format_exc())
    exit(-1)

if SrcTab is None or len(SrcTab)==0:
    print(f"Fatal. Could not load {catchunk} as {catchunk_format}")
    exit(-1)



####################################
# load supporting ELiXer Data
####################################

try:
    elixer_h5 = tables.open_file(elixer_h5_fn)
except:
    print(f"Fatal. Could not load {elixer_h5_fn}",traceback.format_exc())
    exit(-1)

EOTab = Table(elixer_h5.root.ExtractedObjects.read_where("(selected==True) & ((filter_name==b'g')|(filter_name==b'r')|(filter_name==b'f606w'))"))
#note: this will be applied in a loop just a bit further down
#EOTab = EOTab[~eo_sel] #do not restrict here ... this is the general case, not the OII specific case
EOTab.add_index('detectid')

#elixer_h5.close() ... no don't close, need it for reading the lines


################################
# Helpers
################################

def get_source_id_counter():
    """
    get current (to be used) source_id counter and update to the next counter for the next caller

    :return:
    """

    lock = FileLock(source_id_counter_mutex)
    with lock:
        #if it already exists, read in and use the ID, if not, create it and use the default ID
        #either way, incrememnt and write out
        if os.path.exists(source_id_counter_fn):
            current = np.loadtxt(source_id_counter_fn,dtype=np.int64)
            current = np.int64(current) #some weirdness with np.loadtxt, insists on returnedin 0 length array
        else:
            current = global_source_id_counter

        np.savetxt(source_id_counter_fn, [current + 1], fmt='%d')

    return current


def make_selection(flags_cond, flags, selection_type='any'):
    """
    choose selection_type 'any', 'all', 'exact'
    to match any of the flags in the flag_list OR
             all of the flags in the flag_list (the detection must contain all those flags but can also have additional flags)
             exactly and only those flags in the flag_list
    """

    sel_flag = None
    if selection_type == 'any':
        # print("Scanning flags for {total_dets} detections. This may take a while ...")
        sel_flag = np.array([f & flags_cond for f in flags])
    elif selection_type == 'all':
        # print(f"All flag string: {flags_cond:08x}")
        sel_flag = np.array([f & flags_cond == flags_cond for f in flags])
    elif selection_type == 'exact':
        # print(f"Exact flag string: {flags_cond:08x}")
        sel_flag = np.array(flags == flags_cond)
    else:
        print("Invalid selection_type: must be one of: []'any','all','exact'] ")

    return np.array(sel_flag != 0)


def load_source_table(flag_best=True):
    # version = '5.0.1'
    version = '5.0.2'
    catfile = os.path.join(config.hdr_dir['hdr5'], 'catalogs', 'source_catalog_' + version + '.fits')
    print("Loading source table ...")
    st = Table.read(catfile)
    if flag_best:
        st = st[st['flag_best'] == 1]
    print('Source catalog was found at {}'.format(catfile))
    return st


def flatten_list(nested_list):
    def gen_flatten_list(nested_list):
        for item in nested_list:
            if isinstance(item, list) or isinstance(item, np.ndarray):
                yield from gen_flatten_list(item)
            else:
                yield item

    return list(gen_flatten_list(nested_list))


common_rest_lines = np.array([
    G.LyA_rest,
    G.OII_rest,
    G.OIII_4959,
    G.OIII_5007,
    G.OVI_1035,
    G.CIV_1549,
    G.CIII_1909,
    G.CII_2326,
    G.MgII_2799,
    G.Hbeta_4861,
    G.Hgamma_4340,
    G.Hdelta_4101,
    # G.Hepsilon_3970, #others should be there, stop at Hdelta
    # G.Hzeta_3889,
    # G.Heta_3835,
    G.HeI_3888,
    G.HeI_4471,
    G.NV_1241,
    G.SiII_1260,
    G.SiIV_1400,
    G.HeII_1640,
    G.NeIII_3869,
    G.NeIII_3967,
    G.NeV_3347,
    G.NeVI_3427,
    G.NaI_4980,
    G.NaI_5153,
    # absorption
    G.CaII_K_3934,
    G.CaII_H_3968,
])


def line_z_consistency_check(z, found_lines, max_delta_aa=6.0):
    """
    very basic ... is the redshift consistent with any of the found_lines (e.g. assuming a line could be one of 2 dozen common lines)?
    z (float)
    found_lines (array of floats, observed space wavelength in AA)
        Note: this needs to include the ELiXer found lines + the original HETDEX anchor line
    give a little extra room to LyA for velocity offset and broadness
    returns list of rest_wavelengths that are consistent

    Note: could add further checks, along the lines of ELiXer for whether then lines are internally consistent
          that is, do they belong together and are the in the right range of ratios and widths?

    COULD consider changing the input z (if a consistent line is found) to the z based on the line's rest wavelength
          but it may be best to leave as is so that the catalog z_hetdex_src remains strictly True
    """

    # max_delta_aa = 6.0
    # consistent = False

    close_lines = []
    # is the redshift consistent with any?
    for line in found_lines:
        w_rest = line / (1.0 + z)
        w_close = np.isclose(w_rest, common_rest_lines, atol=max_delta_aa)
        if np.any(w_close):
            close_lines.append(common_rest_lines[w_close])
    try:
        # print(f"*** test *** close_lines pre, {close_lines}")
        close_lines = np.array(flatten_list(close_lines))
        # print(f"*** test *** close_lines, post, {close_lines}")
        close_lines = np.array(close_lines).flatten()
    except:
        print(f"Exception in line_z_consistency_check; {traceback.format_exc()}")
    return close_lines


#
# (More) Helpers
#

def evaluate_z_hetdex(row):
    """
      row is the row from the source catalog
      either this is a lone detection OR it is the

      return bool: True (if good), False (if z_hetdex changed)
             float: z_hetdex
             bstr: z_hetdex_src
    """
    try:
        det = row['detectid']
        srcid = row['source_id']
        z = None  # the z we will adjust
        z_src = None
        z_agn = row['z_agn']
        z_diagnose = row['z_diagnose']
        z_elixer = row['z_elixer']
        z_hetdex = row['z_hetdex']
        z_hetdex_src = row['z_hetdex_src']
        flags = row['flags']  # elixer flags
        good = True  # start with this assumption

        # we are normally to trust z_agn as they are set manually
        # 2em are two emission lines that match LyA + CIV and are assigned by Erin ... not a Chenxu manual check AGN
        # is this from Chenxu
        if z_hetdex_src == "liu_agn":
            # print("Chenxu AGN")
            return True, z_hetdex, z_hetdex_src

        # if (z_agn is not None and z_agn > 0) or z_hetdex_src == "liu_agn":
        #     #good as is
        #     if z_hetdex == z_agn:
        #         return True, z_hetdex, z_hetdex_src
        #     else:
        #         return False, z_agn, z_hetdex_src #needs to be updated

        # there are some issues in the data where z_agn does not match z_hetdex but it should

        # just for test?
        # if (z_hetdex != z_elixer) and (z_hetdex != z_diagnose):
        #    print(f"{row['detectid']} unmatched z_hetdex {z_hetdex_src}: {z_hetdex:0.4f} != {z_agn:0.4f} or {z_elixer:0.4f} or {z_diagnose:0.4f}")

        # otherwise, the ASSUMPTION is this came from elixer or diagnose ... anything else would be considered wrong and should be corrected anyway
        if z_diagnose + z_elixer != 0:
            zx = 2. * abs((z_diagnose - z_elixer) / (z_diagnose + z_elixer))  # difference/average
        else:  # above logic prevents large z, so these if sum to zero most both be near zero
            zx = 0

        if zx <= 0.1:  # they are close, so use the elixer value (nominally based on an emission line)
            # use the elixer value
            z = z_elixer
            z_src = 'elixer'
            good = False
        elif row['gmag'] < 23.0:  # normally we use 22 but will give it some room here
            # which one to use?
            # if it is bright, probably z_diagnose
            # or if z_diagnose matches with a line in elixer that could be OII?
            # BUT !!! what about, say a star that really does not have a line? the emission line is bogus?

            lines = elixer_h5.root.SpectraLines.read_where("detectid==det")
            line_waves = np.array(list([row['wave']]) + list(lines['wavelength'])).flatten()

            lines_consistent = line_z_consistency_check(z_diagnose, line_waves, max_delta_aa=6.0)

            if len(lines_consistent) >= 1:
                # keep the z_diagnose but tweak it to the best fit line?
                # print(f"*** test: {row['wave']} {lines_consistent}")
                if row['wave'] >= 3470.0:
                    z_consistent = np.array([w / row['wave'] + 1.0 for w in lines_consistent])
                    # which z_consistent is closest to z_diganose
                    z = z_consistent[np.abs(z_consistent - z_diagnose).argmin()]
                else:  # this is a continuum source, just keep it as is
                    z = z_diagnose
                z_src = 'diagnose'
                good = False
            else:
                # this COULD still be something like a star (continuum source) with a bad line
                # maybe there is imaging? If it is very round and Diagnose says z = 0 (or nearly 0) take it as a star
                if z_diagnose < 0.001:

                    sn = row['sn']
                    chi2 = row['chi2']
                    gmag = row['gmag']
                    eo_rows = EOTab[EOTab['detectid'] == det]
                    for eo_row in eo_rows:
                        if eo_row['filter_name'] == b'f606w':
                            continue  # I don't trust the f606w astrometry we are using to this precision for this purpose
                        # sometimes see odd fits in source extractor, so just check if either set is near circular
                        try:
                            major = eo_row['major']
                            minor = eo_row['minor']
                            e = (major - minor) / major
                        except:
                            e = 999.  # so won't trip the next if
                            major = 999.
                            minor = 999.

                        if (e < 0.1 and major < 6.0) or ((e > 99) and (sn < 6.5) and (chi2 > 1.5) and (gmag < 21.0)):
                            # pretty circular, reasonably pointsource plausible, could still be an elliptical, but may as well assume a star
                            # or bright with weak line
                            # todo: could do better by seeing if can find H & K (common)
                            z = z_diagnose
                            z_src = 'diagnose'
                            good = False
                            # actually, in this case, we'd say NOT a galaxy
                            # sel_eo_dist[i] = False #remove from the list
                            # print(f"remove from list: {det}")
                            break  # no need to check the next one

                    # if sel_eo_dist[i]: #if it is still True, give it the elixer redshift
                    if flags & (G.DETFLAG_BAD_EMISSION_LINE | G.DETFLAG_BAD_FIBERTRACE):
                        z = z_diagnose
                        z_src = 'diagnose'
                        good = False
                        questionable_z_dets.append(det)
                    else:
                        # could still be a problem
                        # bright neighbor as (as a star?) could still confuse it
                        # if the gmag is bright and the redshift if low ( < 0.05) THERE really should be corroborating lines
                        if row['gmag'] < 21 and flags & G.DETFLAG_LARGE_NEIGHBOR:  # this is suspect
                            sel_lines = lines['used'] * lines['continuum'] > 0
                            if np.count_nonzero(sel_lines) < 2:  # really expect OII-3727, both OIII, H-beta, h-gamma
                                # sel_eo_dist[i] = False #remove from the list .. assume this to be outskirts of star
                                z = z_diagnose  # have to give it something so as to not break later
                                z_src = 'diagnose'
                                good = False

                        if z is None:
                            z = z_elixer
                            z_src = 'elixer'
                            good = False
                            questionable_z_dets.append(det)

                else:
                    if flags & (G.DETFLAG_BAD_EMISSION_LINE | G.DETFLAG_BAD_FIBERTRACE):
                        z = z_diagnose
                        z_src = 'diagnose'
                        good = False
                    else:
                        z = z_elixer  # I do question this one ... maybe we should add a flag here? if there are not many we can examine manually
                        z_src = 'elixer'
                        good = False
                        questionable_z_dets.append(det)


        else:  # too faint for Diagnose
            z = z_elixer  # our best bet
            z_src = 'elixer'
            good = False

        if z == z_hetdex:  # it did not change
            return True, z_hetdex, z_hetdex_src
        else:
            return False, z, z_src
    except:
        print(f"[{row['detectid']}]Exception", traceback.format_exc())
        return True, -1.0, None  # exception case, but have to othewise assume z_hetdex is good (that is, do not change it)


def split_cluster(rows, atol=0.005):
    """
    break into 1 or more sub-lists of rows (detections) that should be re-grouped ... new source_ids will be assigned

    this ASSUMES z_hetdex is now correct and it clusters on z_hetdex

    """
    new_cluster = []
    ids = np.zeros(len(rows)).astype(int)
    done = np.full(len(rows), False)

    iterct = 0

    while not all(done):
        uid = np.unique(ids[~done])
        for xid in uid:
            # print(f"({iterct}): ids={ids} uids={uid} xid={xid} done={np.count_nonzero(done)}")
            iterct += 1
            if iterct > 10:
                print(f"WTF? Breaking ... too many iterations. original source_id = {rows['source_id'][0]}")
                done = np.full(len(rows), True)
                break

            # split into 2
            sel = np.array(ids == xid)
            if np.count_nonzero(sel) == 0:
                continue
            rx = rows['z_hetdex'][sel]
            marker = rx[0]
            s1 = np.array([abs(x - marker) < atol for x in rows['z_hetdex']])  # repeats over the whole array
            s0 = ~s1
            # print(f"     {avg} {s1} {s0} {rows['z_hetdex'][sel]}")

            done = done | s1  # if it is every TRUE then it stays TRUE since we are "done" with that entry

            # otherwise update and repeat the loop
            ids[s0 & ~done] += 1  # only update those that are NOT done already

    # build up the clusters by detectid
    uid = np.unique(ids)
    # print(f"Build new lists. uid = {uid}  ids = {ids}")
    for xid in uid:
        sel = np.array(ids == xid)
        new_cluster.append(list(rows['detectid'][sel]))

    return new_cluster


#
# should this be on the whole cluster, not just the current detectid and the seldet for the cluster?
#
def evaluate_cluster(row, row_seldet):
    """
    evaluate if the row (this detection) should be clustered with row_seldet
    where row != row_seldet

    Basically do these two detections belong together?

    shared emission lines
    similar magnitude ? (might not apply if way off the object)
    same EOTab counterpart
    distance between them (esp. given the magnitude)


    """
    keep_cluster = True
    try:

        # first is the row_seldet a good z_hetdex?
        # evaluate each independently
        keep_z_hetdex, new_z_hetdex, new_z_hetdex_src = evaluate_z_hetdex(row)
        seldet_keep_z_hetdex, seldet_new_z_hetdex, seldet_new_z_hetdex_src = evaluate_z_hetdex(row_seldet)

        # if neither changed, assume the grouping is still good?
        if keep_z_hetdex == seldet_keep_z_hetdex == True:
            return True  # keep the cluster

        # one or the other needs to change
        # could be they are both wrong in the same way and still belong together
        # if (keep_z_hetdex == seldet_keep_z_hetdex == False) and \
        #     np.isclose(new_z_hetdex,seldet_new_z_hetdex,atol=0.01):
        #     #todo: they need to change, but maybe still stay clustered
        #     pass

        #######################################
        # check emission lines (type == 1)
        #######################################
        # do they share any emission lines? note that row['wave'] could be 0  for continuum
        det = row['detectid']
        # print("*** det",det)
        lines = elixer_h5.root.SpectraLines.read_where("(detectid==det) & (type==1)")
        add_extra_line = False
        if row['wave'] >= 3470.:
            if any(np.isclose(row['wave'], lines['wavelength'], atol=6.0)):
                line_waves = lines['wavelength']
            else:
                line_waves = list([row['wave']]) + list(lines['wavelength'])
                add_extra_line = True
        else:
            line_waves = lines['wavelength']
            # print(f"*** test ***: row['wave'] == {row['wave']}")

        if add_extra_line:
            line_scores = list([-1]) + list(
                lines['score'])  # use -1 for the passed in line, elixer should have found it too, if it is strong
            # and it would be in the lines list
        else:
            line_scores = lines['score']

            # this is the selected det, so must assume its waves to be "good"
        seldet_det = row_seldet['detectid']
        # print("*** seldet_det",seldet_det)
        seldet_lines = elixer_h5.root.SpectraLines.read_where(
            "(detectid==seldet_det) & (score > 5.0) & (sn > 4.8) & (sigma < 8.0) & (type==1)")
        add_extra_line = False
        if row_seldet['wave'] >= 3470.:
            if any(np.isclose(row_seldet['wave'], seldet_lines['wavelength'], atol=6.0)):
                seldet_line_waves = seldet_lines['wavelength']
            else:
                seldet_line_waves = list([row_seldet['wave']]) + list(seldet_lines['wavelength'])
                add_extra_line = True
        else:
            seldet_line_waves = seldet_lines['wavelength']
            # print(f"*** test ***: row_seldet['wave'] == {row_seldet['wave']}")

        # seldet_line_waves = list([row_seldet['wave']])+list(seldet_lines['wavelength'])
        if add_extra_line:
            seldet_line_scores = list([np.nan]) + list(seldet_lines['score'])  # use nan to ignore the passed in value
        else:
            seldet_line_scores = seldet_lines['score']

        # full cross check, keep it simple though, just checking the wavelength
        # could check on the score as well or the flux maybe (but if moving off source, they could get poor)
        wave_match_ct = 0
        wave_no_match_ct = 0  # what about significant line(s) in det that don't match in selected det?
        swave_match = np.full(len(seldet_line_scores), False)
        swave_no_match_scores = []
        for w, s in zip(line_waves, line_scores):
            if w > 3470.:
                sel_is_close = np.isclose(w, seldet_line_waves, atol=6.0)
                swave_match = swave_match | sel_is_close
                ct = np.count_nonzero(sel_is_close)
                if ct == 0:
                    if s > 8.0:
                        wave_no_match_ct += 1  # significant line in the det that is not found in the selected det
                        swave_no_match_scores.append(s)
                else:
                    wave_match_ct += ct

        seldet_line_scores = np.array(seldet_line_scores)
        if wave_match_ct > 0 and not all(np.isnan(seldet_line_scores[swave_match])):
            # this can be empty???
            # if np.count_nonzero(swave_match) == 0:
            #     print(f"[{row['detectid']}] wtf?? swave_match is 0 length but wave_match_ct == {wave_match_ct}" )
            # if len(seldet_line_scores) == 0:
            #     print(f"[{row['detectid']}] wtf?? seldet_line_scores is 0 length but wave_match_ct == {wave_match_ct}" )
            swave_match_mean = np.nanmean(seldet_line_scores[swave_match])
        else:
            swave_match_mean = 0.
        if len(swave_no_match_scores) > 0 and not all(np.isnan(swave_no_match_scores)):
            swave_no_match_mean = np.nanmean(swave_no_match_scores)
        else:
            swave_no_match_mean = 0.

        if (wave_match_ct == 0 and row['gmag'] > 22.0) or (
                wave_no_match_ct > 0 and swave_no_match_mean > swave_match_mean * 1.33):
            return False

        # clustering ... check that some waves overlap,
        # if imaging, are they on the same counterpart (e.g. check RA, Dec reported  ... either only the selected one or any?)



    except:
        print("Exception", traceback.format_exc())
    return True


# note: this might be a wholly separate loop
def evaluate_counterpart_dist(row):
    """
    """
    out_dict = {
        "counterpart_ra": -999.9,
        "counterpart_dec": -999.9,
        "counterpart_mag": -999.9,
        "counterpart_mag_err": -999.9,
        "counterpart_dist": -999.9,
        "counterpart_catalog_name": "",
        "counterpart_filter_name": "",
        "counterpart_image_depth_mag": -999.9,
    }
    try:

        # EOTab already down selected to be minimum and have g or r band only
        eo_rows = EOTab[EOTab['detectid'] == row['detectid']]

        if len(eo_rows) > 0:
            # use 'g' preferred
            # get closest
            if len(eo_rows) == 1:
                eo_i = 0
            else:
                eo_i = np.argmin(eo_rows['dist_baryctr'])  # these should all be > 0  for this selection

            # print(f"selection 1: index = {eo_i}, dist = {eo_rows['dist_baryctr'][eo_i]}")
            out_dict["counterpart_ra"] = eo_rows['ra'][eo_i]
            out_dict["counterpart_dec"] = eo_rows['dec'][eo_i]
            out_dict["counterpart_mag"] = eo_rows['mag'][eo_i]
            out_dict["counterpart_mag_err"] = eo_rows['mag_err'][eo_i]
            out_dict["counterpart_dist"] = eo_rows['dist_baryctr'][eo_i]
            out_dict["counterpart_catalog_name"] = eo_rows['catalog_name'][eo_i]
            out_dict["counterpart_filter_name"] = eo_rows['filter_name'][eo_i]
            out_dict["counterpart_image_depth_mag"] = eo_rows['image_depth_mag'][eo_i]

        return out_dict
    except:
        print("Exception", traceback.print_exc())
    return out_dict


#####################################
# make sure every entry has a sourceId
# and every soruceID has a selected_det
#####################################

u_srcs, u_cts = np.unique(SrcTab['source_id'], return_counts=True)
sel = u_cts > 1
keep_srcs = u_srcs[sel]
#print(len(keep_srcs))

sel = np.full(len(SrcTab), False)
sel_seldet = np.array(SrcTab['flag_seldet'] != 0)

#print("sub select 1")
SubSrcTab = SrcTab[sel_seldet]

# want the indicies into SubSrcTab that match the source_ids in u_srcs
#print("intersect ...")
_, i1, i2 = np.intersect1d(keep_srcs, SubSrcTab['source_id'], assume_unique=True, return_indices=True)
#print(len(i1), len(i2))

#print("sub select 2")
SubSrcTab = SubSrcTab[i2]

#print(len(SubSrcTab), len(keep_srcs), len(np.intersect1d(keep_srcs, SubSrcTab['source_id'], assume_unique=True)))

no_seldet_sources = np.setdiff1d(keep_srcs, SubSrcTab['source_id'])  # 74 of them

# confirm the 1:1 match ... huh they are not 1:1  ...close but not quite
# 436984 437058 436984 ... there are 74 more sources in keep_srcs than are found in SubSrcTab
# why? these 74 do NOT have any flag_seldet in the source_catalog ... no detectid is set as the seldet
# So, what do we want to do??? well, lets pick one to be the seldet for now and re-insert them
if len(no_seldet_sources) > 0:
    print(f"[{catchunk}] re-assigning a selected_det for those that do not have one ...")
    no_seldet_sources = np.setdiff1d(keep_srcs, SubSrcTab['source_id'])
    no_seldet_set_det = []
    for s in no_seldet_sources:
        rows = SrcTab[SrcTab['source_id'] == s]

        # need to pick a seldet
        # of those that are emission lines, pick the one with the highest line flux
        # if all are continuum, pick the largest gmag
        s_line = np.array(rows['det_type'] == 'line')
        if np.count_nonzero(s_line) > 0:
            ix = np.argmax(rows['flux'][s_line])
            d = rows['detectid'][s_line][ix]
        else:
            ix = np.argmax(rows['gmag'])
            d = rows['detectid'][ix]

        no_seldet_set_det.append(d)
        sel = np.array(SrcTab['detectid'] == d)
        SrcTab['flag_seldet'][sel] = 1
        SrcTab['selected_det'][sel] = True

    #print(f"[{catchunk}] re-select SubSrcTab ...")

    sel = np.full(len(SrcTab), False)
    sel_seldet = np.array(SrcTab['flag_seldet'] != 0)

    #print(f"[{catchunk}] sub select 1")
    SubSrcTab = SrcTab[sel_seldet]

    # want the indicies into SubSrcTab that match the source_ids in u_srcs
    #print(f"[{catchunk}] intersect ...")
    _, i1, i2 = np.intersect1d(keep_srcs, SubSrcTab['source_id'], assume_unique=True, return_indices=True)
    #print(len(i1), len(i2))

    #print(f"[{catchunk}] sub select 2")
    SubSrcTab = SubSrcTab[i2]

    #print(len(SubSrcTab), len(keep_srcs), len(np.intersect1d(keep_srcs, SubSrcTab['source_id'], assume_unique=True)))
    print(f"[{catchunk}] SrcTab (as selected) has been modified s|t every source_id now has exactly one detectid that is the selected_det")

print(f"[{catchunk}] Now SubSrcTab has one entry for every source_id and that entry is the detectid that is the selected_det")


#############################
# main loop
#############################

#
# a check of z_hetdex both in terms of detectid (single detetctid for a source_id) and as a result of clustering
#



questionable_z_dets = []  # list of detectids that we really are not certain about the redshift
changed_z_dets = []  # dets where z_hetdex changed
changed_counterpart_dets = []
# this can ONLY turn "True" selection into False
# note: if you want this independent of the other selections, then
#      need to open this up to all SrcTab: i.e. define a new sel_xxx = np.full(len(SrcTab),True) and iterate over that
# for i,row in enumerate(tqdm(SrcTab[sel])):

# choose a specifc one to check
sel_row = None  # SrcTab['detectid'] == 4018105236

# z to consider
# z_agn
# z_elixer
# zspec
# z_diagnose
# z_hetdex
# z_hetdex_conf ??
# z_hetdex_src  ... who got it, one of (byte strings):   0.0, 2em, diagnose, elixer, liu_agn

dets_with_errors = []
dets_with_errors_reason = []



source_id_mismatches = []

# sdet = SrcTab['detectid'] == check_det

# for i,row in enumerate(tqdm(SrcTab[sdet])):
for i, row in tqdm(enumerate(SrcTab),total= len(SrcTab),disable=not SHOW_TQDM):

    if check_counterpart_dist:
        cp_dict = evaluate_counterpart_dist(row)
        if (cp_dict is not None) and \
                (row['counterpart_dist'] != cp_dict['counterpart_dist'] or row['counterpart_mag'] != cp_dict[
                    'counterpart_mag']):
            sel_src_det = SrcTab['detectid'] == row['detectid']
            SrcTab['counterpart_ra'][sel_src_det] = cp_dict['counterpart_ra']
            SrcTab['counterpart_dec'][sel_src_det] = cp_dict['counterpart_dec']
            SrcTab['counterpart_mag'][sel_src_det] = cp_dict['counterpart_mag']
            SrcTab['counterpart_mag_err'][sel_src_det] = cp_dict['counterpart_mag_err']
            SrcTab['counterpart_dist'][sel_src_det] = cp_dict['counterpart_dist']
            SrcTab['counterpart_catalog_name'][sel_src_det] = cp_dict['counterpart_catalog_name']
            SrcTab['counterpart_filter_name'][sel_src_det] = cp_dict['counterpart_filter_name']
            SrcTab['counterpart_image_depth_mag'][sel_src_det] = cp_dict['counterpart_image_depth_mag']
            changed_counterpart_dets.append(row['detectid'])

    if not check_z_hetdex:
        continue

    det = row['detectid']

    # print(f"\n *** test *** starting {det} index {i}")

    z_elixer = row['z_elixer']
    z_diagnose = row['z_diagnose']
    keep_z_hetdex = True  # assume z_hetdex is correct to start

    # if we want to change it, need to check to see if this detectid is in a cluster and if so
    #   do we need to remove it from the cluster
    #     and if we do remove it, does it belong with any other detectids in a different cluster
    #        (so maybe if we remove it, check the others in the same original cluster for a wavelength match and move them all?)

    if z_diagnose and z_elixer <= -1.0:
        # nothing to fall back on, have to accept it as it is
        print(f"{det} {z_diagnose} {z_elixer} ; cannot adjust, no reference")
        continue

    # passed the initial check, we do have something to refence, so load the rest of the data we want
    srcid = row['source_id']
    flags = row['flags']  # elixer flags
    z = None  # the z we will adjust
    z_agn = row['z_agn']
    z_hetdex = row['z_hetdex']
    z_hetdex_src = row['z_hetdex_src']

    # first check if z_hetdex agrees with z_elixer and z_diagnose and z_agn (or they are < -1, not to be used)
    # then we assume good and move on
    # that is, every contributor to z agrees or is not set
    if ((z_agn < 0) or (np.isclose(z_agn, z_hetdex, atol=atol_strong))) and \
            ((z_diagnose < 0) or (np.isclose(z_diagnose, z_hetdex, atol=atol_weak))) and \
            ((z_elixer < 0) or (np.isclose(z_elixer, z_hetdex, atol=atol_strong))):
        # print(f"{det} PASS {z_hetdex:0.4f} {z_agn:0.4f} {z_diagnose:0.4f} {z_elixer:0.4f}")
        continue
        # else:
    #    #just for testing
    #    print(f"{det} FAIL {z_hetdex:0.4f} {z_agn:0.4f} {z_diagnose:0.4f} {z_elixer:0.4f}")

    # ok. there is a mismatch of sorts
    # easiest case is that we have the z_diagnose redshift and it is just not as precise as z_elixer
    # BUT could be a bad clustering, or a method disagreement
    #  z_hetdex_src values:  0.0, 2em, diagnose , elixer , liu_agn
    # simple rules ... if NOT clustered (e.g. if this source_id has only one detections)
    # if z_agn is set or the z_hetdex_src is "liu_agn", it wins. We keep z_hetdex and move on (noteL 2em is an auto-set by Erin if LyA + CIV)
    # if z_hetdex

    # print(f"iter: {i}")

    sel_src = np.array(SrcTab['source_id'] == srcid)
    src_ct = np.count_nonzero(sel_src)
    det_seldet = None  #by default only the one, so don't need to set this(iu\f needed, is set later)
    row_z_hetdex = None
    if src_ct == 1:  # this is not clustered, simple case
        # BUT z_agn might not be set ... could be -1
        if z_hetdex_src == 'liu_agn':  # z_agn is not necessarily from Chenxu
            # if (z_hetdex_src == 'liu_agn' and (np.isclose(z_agn,z_hetdex,atol=atol_strong)):
            # print(f"{det} Chenxu AGN")
            continue  # this is good as is

        # otherwise there is disagrement and we need to check
        z = None
        row_z_hetdex = None

    elif src_ct == 0:  # this is an error case, should not happen
        print(f"[{catchunk}] {det} WFT? has not source_id??")
        dets_with_errors.append(det)
        dets_with_errors_reason.append("no source_id")
        continue  # have to just ignore it and move on to the next
    else:
        # this is clustered and the mismatch COULD be due to a bad clustering
        # either it does not belong in the cluster OR the cluster has a bad redshift
        # print(f"{det} clustered")

        # for reference which one is the selected det
        # this is a complicated subselection of a subselection BE CAREFUL
        which_is_seldet = SrcTab['flag_seldet'][sel_src] != 0
        if np.count_nonzero(which_is_seldet) != 1:
            # this is a problem ... most likely due to early cuts, the one with the flag_seldet was excluded
            # this PROBABLY means that it is a poorly chosen seldet for this clustering
            print(
                f"[{catchunk}] [{det}] Error. Source_id {srcid} has no flag_seldet. Claims {src_ct} members. selection = {np.count_nonzero(which_is_seldet)}")

            # alternate check
            which_is_seldet = SrcTab['selected_det'][sel_src] == True
            if np.count_nonzero(which_is_seldet) != 1:
                print(
                    f"[{catchunk}] [{det}] Error. Alternate check: Source_id {srcid} has no flag_seldet. Claims {src_ct} members. selection = {np.count_nonzero(which_is_seldet)}")
                continue

            # if z_hetdex_src == 'liu_agn':  # z_agn is not necessarily from Chenxu
            #     continue  # this is good as is
            #
            # # otherwise there is disagrement and we need to check
            # z = None
            # row_z_hetdex = None
        else:
            det_seldet = SrcTab[sel_src][which_is_seldet]['detectid']
            # print(f"*** test: det_seldet = {det_seldet}")

        # which det in the cluster is the source of this redshift? note that it need not be the one that is the "flag_seldet"
        # what if there is more than one?
        s1 = np.array(SrcTab['z_hetdex'][sel_src] == z_hetdex)
        s1_ct = np.count_nonzero(s1)
        if s1_ct == 0:  # this is an error case, no one has this z ??
            print(f"{det} [{srcid}] Error! none have matching z_hetdex = {z_hetdex}")
            # what do we want to do ?
            row_z_hetdex = None
            dets_with_errors.append(det)
            dets_with_errors_reason.append("none match z_hetdex")
            # we do NOT know the origin of this z_hetdex ... so just re-evaluate it?
            # could it be that it was set manually?
            print("todo: need to rebuild the source_id ...")
        elif s1_ct > 1:  # which one to use?
            # again, it might be that none of these are the flag_seldet
            # print(f"{det} {s1_ct} matches to z_hetdex in source_id")
            try:
                # print(f"*** test *** {len(sel_src)} {np.count_nonzero(sel_src)} {len(s1)} {np.count_nonzero(s1)} {det_seldet}")
                s2 = np.array(SrcTab['detectid'][sel_src][s1] == det_seldet)
                row_z_hetdex = SrcTab[sel_src][s1][s2]
                if len(row_z_hetdex) != 1:  # this is an error
                    print(f"{det} [{srcid}] ERROR! multiple flag_seldet")
                    # do we give up here?
                    dets_with_errors.append(det)
                    dets_with_errors_reason.append("multiple flag_seldet")
                    continue
                else:
                    row_z_hetdex = row_z_hetdex[0]  # take out of the list
            except:
                print(f"Exception: det = {det}, det_seldet == {det_seldet}: {traceback.format_exc()}")
        else:  # ct is exactly 1, use this one
            row_z_hetdex = SrcTab[sel_src][s1][0]
            # print(f"{det} exactly 1 match to z_hetdex in source_id")

    # okay so we have the detectid we want to investigate and (possibly) the row (row_z_hetdex) that is the source of the z_hetdex value
    # row is THIS row, row_z_hetdex is the source of the z_hetdex IF this detection is clustered with other detections
    #    and in that case it could be that row == row_z_hetdex or it could be another detection

    # so options:
    #  1. this is a single, lone detection OR could not find a z_hetdex match in cluster. Either way, treat as an individual
    #     row is populated, row_z_hetdex is None, src_ct = 1
    #  2. this is in a cluster but could not find a z_hetdex match in the cluster.
    #     row is populated, row_z_hetdex is None, src_ct > 1
    #  3. this is in a cluster and we found a single z_hetdex match (or selected the flag_seldet one) in the cluster.
    #     row is populated, row_z_hetdex is populated, src_ct > 1
    #  4. it already looks correct, we don't get this far and we just continue (earlier)
    #
    #  5. errors that we will have to handle later ... we don't get this far. The error list is populated and continue is triggered

    # first if there is a row_z_hetdex, should we be clustered with it?
    #  if row != row_z_hetdex (and row_z_hetdex is not None)
    #    evaluate:  should we be clustered with it?
    #    if NO, uncluster
    #      then re-evaluate the z for THIS row (needs to be its own function)
    #    else
    #      keep the z_hetdex and the cluster and move one (*note: later when we hit the row == row_z_hetdex, it will may get updated)
    #  else: # this row IS the z_hetdex_row
    #    re-evaluate the z for THIS row
    #    did the z_hetdex change?
    #    if Yes
    #      find all other detections on this source_id and update their z_hetdex
    #    else:
    #      do nothing

    if (row_z_hetdex is not None) and (row['detectid'] != row_z_hetdex['detectid']):
        # this is a clustering and THIS detection is NOT the seldet
        # print("todo: evaluate_cluster()")
        keep_cluster = evaluate_cluster(row, row_z_hetdex)

        if keep_cluster:  # fine as is
            # print(f"{det} KEEP cluster")
            pass
        else:
            # need to split up?
            # print(f"{det} need to uncluster")

            # mark this Source_id to be re-clustered (keep a list of source_ids to be re-worked)
            #  but don't change it yet, will do that later as a separate loop as it is not
            #  sufficient to break this detectid out of the source, we may need to split the old source up into
            #  multiple new sources

            source_id_mismatches.append(srcid)

            # if srcid not in source_id_mismatches:
            #    #or could blindly add then do unique at the end
            #    source_id_mismatches.append(srcid)

            keep_z_hetdex, new_z_hetdex, new_z_hetdex_src = evaluate_z_hetdex(row)
            if not keep_z_hetdex:  # it changed
                # update the source_id (might just be a single detection or may be multiples)
                # remeber this is either a single OR it is the seldet for the source
                SrcTab['z_hetdex'][sel_src] = new_z_hetdex
                SrcTab['z_hetdex_src'][sel_src] = new_z_hetdex_src
                changed_z_dets.append(det)

    else:  # this is a single detection OR it is the cluster's seldet
        # print("Calling evaluate")
        keep_z_hetdex, new_z_hetdex, new_z_hetdex_src = evaluate_z_hetdex(row)
        if not keep_z_hetdex:  # it changed
            # update the source_id (might just be a single detection or may be multiples)
            # remeber this is either a single OR it is the seldet for the source

            update_sel_det = np.array(SrcTab['detectid'] == det)
            SrcTab['z_hetdex'][update_sel_det] = new_z_hetdex
            SrcTab['z_hetdex_src'][update_sel_det] = new_z_hetdex_src
            changed_z_dets.append(det)

            # AND update the others in the source only IF they match the old z_hetdex
            # now it could be that we have not gotten to some of these yet and they will end up being un-clustered,
            # but if we hit them later, they will still have their redshift evaluated and can get updated z_hetdex then
            #  and then when we run the re-clustering, they can be broken out
            update_sel_src = sel_src & np.array(SrcTab['z_hetdex'] == z_hetdex)
            # print(f"{det} changing {np.count_nonzero(update_sel_src)} others in source. From {z_hetdex:0.4f} to {new_z_hetdex:0.4f}")
            SrcTab['z_hetdex'][update_sel_src] = new_z_hetdex
            SrcTab['z_hetdex_src'][update_sel_src] = new_z_hetdex_src

        # else: it did not change so do nothing

source_id_mismatches = np.unique(source_id_mismatches)
# changed_counterpart_dets
print(f"[{catchunk}] Changed continuum counterpart: {len(changed_counterpart_dets)}")
print(f"[{catchunk}] Changed z assignments: {len(changed_z_dets)}")
print(f"[{catchunk}] Questionable z assignments: {len(questionable_z_dets)}")
# print(f"Changed: {changed_z_dets}")
print(f"[{catchunk}] Source IDs to update ({len(source_id_mismatches)}): {source_id_mismatches}")

####################################
#update the source_ids
###################################
# now, go through those mismatches and assign new source_ids
try:
    SrcTab.remove_indices('source_id')
except:
    pass

for i in tqdm(range(len(source_id_mismatches)),disable=not SHOW_TQDM):
    src = source_id_mismatches[i]
    rows = SrcTab[SrcTab['source_id'] == src]
    new_clusters = split_cluster(rows)
    for clu in new_clusters:
        # each one is a list
        global_source_id_counter += 1
        r0 = None
        comp_col = 'sn'

        for det in clu:
            SrcTab['source_id'][SrcTab['detectid'] == det] = get_source_id_counter()
            if r0 is None:
                r0 = SrcTab[SrcTab['detectid'] == det]
                if r0['z_hetdex'] == r0['z_elixer']:
                    comp_col = 'sn'
                else:
                    comp_col = 'continuum'
            else:
                # compare emission line snr? what if there is no line? maybe mag_g_wide (should always exist)
                r1 = SrcTab[SrcTab['detectid'] == det]
                if r0[comp_col] < r1[comp_col]:
                    r0 = r1
        if r0 is not None:
            for det in clu:
                if r0['detectid'] == det:
                    SrcTab['selected_det'][SrcTab['detectid'] == det] = True
                    SrcTab['flag_seldet'][SrcTab['detectid'] == det] = 1
                else:
                    SrcTab['selected_det'][SrcTab['detectid'] == det] = False
                    SrcTab['flag_seldet'][SrcTab['detectid'] == det] = 0

# re-do the index
SrcTab.add_index('source_id')

