"""

make a table of statistics for finding bad fibers from the cmbf.fits files


todd:

"""


import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import os.path as op
import glob
from tqdm import tqdm
import traceback

# lets build a table of fibers over time
# need to record the full fiber information* (specid, slot, ifuid, amp, fiber number, date (or month))
#   *initially might only have the slot
# and track the same stats info and status:
#  mean, median, stddev, ct_0, ct_bad(low), ct_iffy(low), ct_high, status, reason
#
# each row would be one fiber (out of ~35K) for that month, so we'd end up with
#  something like 4 million records by the end of 2026 ... that should actually be okay
#

# can build one from each dictionary, then vstack as needed

# cmbf fiber table
# fT = Table(dtype=[('specid', int), ('ifuslot', int),('ifuid', int),('amp','U2'),('yyyymm',int),('cmbf','U16'),
#                   ('mean',float),('median',float),('stddev',float),
#                   ('ct_0',int),('ct_bad_low',int),('ct_bad_medlow',int),('ct_bad_high',int),
#                   ('status',int),('reason',int)])

# fiberid like 008_093_054_LU_011

def new_fT(stats_dict=None):
    fT = Table(dtype=[('yyyymm', int), ('cmbf', 'U16'), ('fiberid', 'U18'),
                      ('specid', int), ('ifuslot', int), ('ifuid', int), ('amp', 'U2'),
                      ('fibnum', int), ('mean', float), ('median', float), ('stddev', float),
                      ('ct_0', int), ('ct_bad_low', int), ('ct_bad_medlow', int), ('ct_bad_high', int),
                      ('status', int), ('reason', int)])

    if stats_dict is not None:

        try:
            num_fibers = len(stats_dict['fibnum'])

            for _ in range(num_fibers):  # should be 112
                fT.add_row(None)

            # all entriesfor this are the same

            if stats_dict['multiframe'] is not None:
                toks = stats_dict['multiframe'].split("_")
                specid = int(toks[1])
                ifuslot = int(toks[2])
                ifuid = int(toks[3])
            else:
                specid = 0  # don't know this yet
                ifuslot = int(stats_dict['ifu_amp'][1:4])
                ifuid = 0  # don't know this yet

            amp = stats_dict['ifu_amp'][5:7]
            mf_base = f"{str(specid).zfill(3)}_{str(ifuslot).zfill(3)}_{str(ifuid).zfill(3)}"
            cmbf = f"{stats_dict['ifu_amp']}cmbf.fits"

            fT['yyyymm'] = np.full(num_fibers, int(stats_dict['month']))
            fT['cmbf'] = np.full(num_fibers, cmbf)
            fT['specid'] = np.full(num_fibers, specid)
            fT['ifuslot'] = np.full(num_fibers, ifuslot)
            fT['ifuid'] = np.full(num_fibers, ifuid)
            fT['amp'] = np.full(num_fibers, amp)

            fiberid = [f"{mf_base}_{amp}_{str(fnum).zfill(3)}" for fnum in stats_dict['fibnum']]

            # now the per fiber bits
            fT['fiberid'] = fiberid
            fT['fibnum'] = stats_dict['fibnum']
            fT['median'] = stats_dict['median']
            fT['stddev'] = stats_dict['stddev']
            fT['ct_0'] = stats_dict['ct_0']
            fT['ct_bad_low'] = stats_dict['ct_bad']
            fT['ct_bad_medlow'] = stats_dict['ct_iffy']
            fT['ct_bad_high'] = stats_dict['ct_high']
            fT['status'] = stats_dict['status']
            fT['reason'] = stats_dict['reason']
        except:

            try:

                if 'fibnum' in stats_dict.keys():
                    if stats_dict['fibnum'] is not None:
                        num_fibers = len(stats_dict['fibnum'])
                    else:
                        num_fibers = 0
                else:
                    num_fibers = 0

                if num_fibers > 0:
                    for _ in range(num_fibers):  # should be 112
                        fT.add_row(None)
                else: #just add one row
                    fT.add_row(None)
                    num_fibers = 1 #just so an insert happens

                if stats_dict['multiframe'] is not None:
                    toks = stats_dict['multiframe'].split("_")
                    specid = int(toks[1])
                    ifuslot = int(toks[2])
                    ifuid = int(toks[3])
                else:
                    specid = 0  # don't know this yet
                    ifuslot = int(stats_dict['ifu_amp'][1:4])
                    ifuid = 0  # don't know this yet

                amp = stats_dict['ifu_amp'][5:7]
                mf_base = f"{str(specid).zfill(3)}_{str(ifuslot).zfill(3)}_{str(ifuid).zfill(3)}"
                cmbf = f"{stats_dict['ifu_amp']}cmbf.fits"

                fT['yyyymm'] = np.full(num_fibers, int(stats_dict['month']))
                fT['cmbf'] = np.full(num_fibers, cmbf)
                fT['specid'] = np.full(num_fibers, specid)
                fT['ifuslot'] = np.full(num_fibers, ifuslot)
                fT['ifuid'] = np.full(num_fibers, ifuid)
                fT['amp'] = np.full(num_fibers, amp)

                fT['status'] = -999
                fT['reason'] = 0xffffffff

            except:
                print(f"Fail. Could not insert dictionary into table {stats_dict}\n{traceback.format_exc()}")

    # could add other indicies, but for now, these two are sufficient
    fT.add_index('ifuslot')
    fT.add_index('yyyymm')

    return fT



def get_nearest_past_date(date_array, date):
    """
    """
    if len(str(date)) == 6:
        date = int(str(date) + "15")  # we'll just assume mid month for now (don't need to worry about 28-31)

    if type(date_array) == list:
        date_array = np.array(date_array)

    idx = (np.abs(date_array - date)).argmin()

    if date_array[idx] <= date:
        pass  # this is good, it is either the same or the nearest is less than the current date
    else:  # nearest is > current date, so back up one
        idx = idx - 1

    if idx < 0:
        return None
    else:
        return date_array[idx]


def get_ifu_mapping_table():
    """
    """
    ifu_tab_path = "/work/03261/polonius/hetdex/single_shot/ifu_mf_map_full_hetdex.dat"
    ifu_tab = None
    try:
        ifu_tab = Table.read(ifu_tab_path, format="ascii")
    except:
        print(f"Could not load ifu mapping table: {ifu_tab_path}")

    return ifu_tab


def ifu_date_to_multiframe(map_tab, date, ifu):
    """
    for a given date and IFU SLot, return the multiframe (SPECID, IFUSlot, IFUID)
    date is YYYYMMDD

    columns:
    yyyymmdd yyyymm specid ifuslot ifuid
    (int)    (int)  (int)   (int)    (int)

    """

    dlist = []
    d = {"specid": None,
         "ifuslot": None,
         "ifuid": None,
         "ifu_str": None,
         "src_date": None
         }

    yyyymm = None
    yyyymmdd = None

    if map_tab is None:
        map_tab = get_ifu_mapping_table()

    if map_tab is None:
        print(f"Could not map {date} {ifu}. No mapping table.")
        return [], 0

    try:
        ifu_target = int(ifu)

        # date might be yyyymm or yyyymmdd
        if len(str(date)) == 6:  # yyyymm
            yyyymm = int(date)
            sel = np.array(map_tab["yyyymm"] == yyyymm) * np.array(map_tab["ifuslot"] == ifu_target)

        elif len(str(date)) == 8:
            yyyymmdd = int(date)
            sel = np.array(map_tab["yyyymmdd"] == yyyymmdd) * np.array(map_tab["ifuslot"] == ifu_target)
        else:
            print(f"Could not map {date} {ifu}. Bad date.")
            return [], 0

        if np.count_nonzero(sel) == 1:  # perfect, easy
            specid = str(map_tab['specid'][sel][0]).zfill(3)
            ifuslot = str(map_tab['ifuslot'][sel][0]).zfill(3)
            ifuid = str(map_tab['ifuid'][sel][0]).zfill(3)
            ifu_str = f"{specid}_{ifuslot}_{ifuid}"

            d = {"specid": specid,
                 "ifuslot": ifuslot,
                 "ifuid": ifuid,
                 "ifu_str": ifu_str,
                 "src_date": int(map_tab['yyyymmdd'][sel][0])
                 }
            dlist.append(d)

        elif np.count_nonzero(sel) > 1:  # more than one match, can happen for a month
            if yyyymmdd is not None:  # this should not be ... if an exact full date is specified
                print(
                    f"Could not map {date} {ifu}. Unexpected, multiple matches ({np.count_nonzero(sel)}) found for date.")
            else:  # this was just a month check

                strmatch = []

                for row in map_tab[sel]:
                    specid = str(row['specid']).zfill(3)
                    ifuslot = str(row['ifuslot']).zfill(3)
                    ifuid = str(row['ifuid']).zfill(3)
                    ifu_str = f"{specid}_{ifuslot}_{ifuid}"

                    if ifu_str not in strmatch:  # don't have this one yet
                        strmatch.append(ifu_str)
                        d = {"specid": specid,
                             "ifuslot": ifuslot,
                             "ifuid": ifuid,
                             "ifu_str": ifu_str,
                             "src_date": int(row['yyyymmdd'])
                             }
                        dlist.append(d)


        else:  # problem none found
            # find nearest previous match
            if yyyymmdd is not None:
                yyyymmdd = get_nearest_past_date(map_tab['yyyymmdd'], yyyymmdd)
            else:
                yyyymmdd = get_nearest_past_date(map_tab['yyyymmdd'], yyyymm)

            #print(f"nearest previous date: {yyyymmdd}")
            if yyyymmdd is not None:
                sel = np.array(map_tab["yyyymmdd"] == yyyymmdd) * np.array(map_tab["ifuslot"] == ifu_target)
                # print(f"DEBUG2: {yyyymmdd} {ifu}")
                if np.count_nonzero(sel) == 1:  # perfect, easy
                    specid = str(map_tab['specid'][sel][0]).zfill(3)
                    ifuslot = str(map_tab['ifuslot'][sel][0]).zfill(3)
                    ifuid = str(map_tab['ifuid'][sel][0]).zfill(3)
                    ifu_str = f"{specid}_{ifuslot}_{ifuid}"

                    d = {"specid": specid,
                         "ifuslot": ifuslot,
                         "ifuid": ifuid,
                         "ifu_str": ifu_str,
                         "src_date": int(map_tab['yyyymmdd'][sel][0])
                         }
                    dlist.append(d)
                else:  # unexpected
                    # print(f"xx {yyyymmdd} {ifu}")
                    print(
                        f"Could not map {date} {ifu}. Unexpected ({np.count_nonzero(sel)}) matches found for date.")
            else:
                print(f"Could not map {date} {ifu}. No match found.")

        return dlist, len(dlist)



    except:
        print(f"Could not map {date} {ifu}")


def check_for_bad_fibers(basedir, yyyymm, ifu, amp, fiber_list=[]):
    """
    """
    IFFY_FIBER_CMBF_THRESH = 0.4
    BAD_FIBER_CMBF_THRESH = 0.2
    out_log = []
    try:
        fn = op.join(basedir, f"{yyyymm}/i{str(ifu).zfill(3)}a{amp}cmbf.fits")
        with fits.open(fn) as hdu:
            per_fiber_avg = np.nanmedian(hdu[0].data, axis=1)

            if len(fiber_list) > 0:
                # remember fiber list is the fiber number 1 to 112 and we want the index 0 - 111
                idx = np.array([f - 1 for f in fiber_list])
            else:
                idx = np.where(per_fiber_avg <= IFFY_FIBER_CMBF_THRESH)[0]

            for i in idx:
                if per_fiber_avg[i] <= BAD_FIBER_CMBF_THRESH:
                    logstr = f" BAD FIBER: {yyyymm}-i{ifu}{amp}: fib# {str(i + 1).zfill(3)} = {per_fiber_avg[i]:0.4f}"
                elif per_fiber_avg[i] <= IFFY_FIBER_CMBF_THRESH:
                    logstr = f"IFFY FIBER: {yyyymm}-i{ifu}{amp}: fib# {str(i + 1).zfill(3)} = {per_fiber_avg[i]:0.4f}"
                else:
                    logstr = f"OKAY FIBER: {yyyymm}-i{ifu}{amp}: fib# {str(i + 1).zfill(3)} = {per_fiber_avg[i]:0.4f}"

                print(logstr)
                out_log.append(logstr)
    except:
        print(f"Excpetion in check_for_bad_fibers(). {traceback.format_exc()}")

    return out_log


def get_all_bad_fibers_in_cmbf(cmbf_file):
    """
    """
    IFFY_FIBER_CMBF_THRESH = 0.4
    BAD_FIBER_CMBF_THRESH = 0.2
    HIGH_FIBER_CMBF_THRESH = 1.8
    list_ifu_amp = []  # only one, but repeat for alignment
    list_month = []  # only one, but repeat for alignment
    list_fibers = []
    list_values = []
    try:
        # "/scratch/projects/hetdex/lib_calib/yyyymm/iXXXaYYYcmbf.fits"
        toks = cmbf_file.split("/")
        month = toks[-2]
        ifu_amp = toks[-1][:7]

        with fits.open(cmbf_file) as hdu:

            if hdu[0].data is None or np.size(hdu[0].data) < (112 * 1032):
                print(f"Bad data payload, skipping {cmbf_file}")
            else:
                per_fiber_avg = np.nanmedian(hdu[0].data, axis=1)

                idx = np.where(per_fiber_avg <= IFFY_FIBER_CMBF_THRESH)[0]

                for i in idx:
                    # 0 is already flagged as bad
                    if per_fiber_avg[i] <= IFFY_FIBER_CMBF_THRESH and per_fiber_avg[i] != 0.0:
                        list_fibers.append(i + 1)
                        list_values.append(per_fiber_avg[i])
                        list_month.append(month)
                        list_ifu_amp.append(ifu_amp)
    except:
        print(f"Exception in get_all_bad_fibers_in_cmbf(). {cmbf_file}\n{traceback.format_exc()}")

    return list_month, list_ifu_amp, list_fibers, list_values


def get_all_high_fibers_in_cmbf(cmbf_file):
    """
    """
    # IFFY_FIBER_CMBF_THRESH = 0.4
    # BAD_FIBER_CMBF_THRESH = 0.2
    HIGH_FIBER_CMBF_THRESH = 1.8
    list_ifu_amp = []  # only one, but repeat for alignment
    list_month = []  # only one, but repeat for alignment
    list_fibers = []
    list_values = []
    try:
        # "/scratch/projects/hetdex/lib_calib/yyyymm/iXXXaYYYcmbf.fits"
        toks = cmbf_file.split("/")
        month = toks[-2]
        ifu_amp = toks[-1][:7]

        with fits.open(cmbf_file) as hdu:

            if hdu[0].data is None or np.size(hdu[0].data) < (112 * 1032):
                print(f"Bad data payload, skipping {cmbf_file}")
            else:
                per_fiber_avg = np.nanmedian(hdu[0].data, axis=1)

                idx = np.where(per_fiber_avg >= HIGH_FIBER_CMBF_THRESH)[0]

                for i in idx:
                    # 0 is already flagged as bad
                    # if per_fiber_avg[i] <= IFFY_FIBER_CMBF_THRESH and per_fiber_avg[i] != 0.0:
                    if True:
                        list_fibers.append(i + 1)
                        list_values.append(per_fiber_avg[i])
                        list_month.append(month)
                        list_ifu_amp.append(ifu_amp)
    except:
        print(f"Exception in get_all_bad_fibers_in_cmbf(). {cmbf_file}\n{traceback.format_exc()}")

    return list_month, list_ifu_amp, list_fibers, list_values


def get_cmbf_fiber_avg(cmbf_file, fibernum):
    """
    """
    try:
        avg = None
        with fits.open(cmbf_file) as hdu:
            per_fiber_avg = np.nanmedian(hdu[0].data, axis=1)

        avg = per_fiber_avg[fibernum - 1]

    except:
        print(f"Exception in get_cmbf_fiber_avg(). {cmbf_file}\n{traceback.format_exc()}")

    return avg


def get_cmbf_fiber_stats(cmbf_file, fibernum,map_tab=None):
    """
    """
    IFFY_FIBER_CMBF_THRESH = 0.4
    BAD_FIBER_CMBF_THRESH = 0.2
    HIGH_FIBER_CMBF_THRESH = 1.8

    d = {}

    if not 1 <= fibernum <= 112:
        d['log'] = "invalid fiber number"
        return d

    if not op.exists(cmbf_file):
        d['log'] = "file does not exist"
        return d

    try:
        with fits.open(cmbf_file) as hdu:
            fiber = hdu[0].data[fibernum - 1, :]

        toks = cmbf_file.split("/")
        month = toks[-2]
        ifu_amp = toks[-1][:7]

        d['month'] = month
        d['ifu_amp'] = ifu_amp

        mflist, mfct = ifu_date_to_multiframe(map_tab, month, toks[-1][1:4])

        if mfct > 0:
            # todo: for now just use the first one, if there is more ... or should we use the last??
            # or should we use the one with the longest time for the month
            d['multiframe'] = f"multi_{mflist[0]['ifu_str']}_{ifu_amp[-2:]}"
        else:
            d['multiframe'] = None

        d['fibnum'] = fibernum  # (fibernum-1) % 112 + 1 #use this IF we allow any value and just modulo 112

        d['len'] = len(fiber)
        d['mean'] = np.nanmean(fiber)
        d['median'] = np.nanmedian(fiber)
        d['stddev'] = np.nanstd(fiber)

        d['ct_0'] = np.count_nonzero(fiber == 0)
        d['ct_bad'] = np.count_nonzero(fiber <= BAD_FIBER_CMBF_THRESH)
        d['ct_iffy'] = np.count_nonzero(fiber <= IFFY_FIBER_CMBF_THRESH) - d['ct_bad']
        d['ct_high'] = np.count_nonzero(fiber >= HIGH_FIBER_CMBF_THRESH)

        d['status'] = -999  # unset (unchecked)
        d['reason'] = 0  # unset
        d['log'] = "ok"
    except:
        d['log'] = traceback.format_exc()
        print(f"Exception in get_cmbf_fiber_avg(). {cmbf_file}")

    return d


def get_cmbf_all_fiber_stats(cmbf_file, map_tab=None):
    """
    """
    IFFY_FIBER_CMBF_THRESH = 0.4
    BAD_FIBER_CMBF_THRESH = 0.2
    HIGH_FIBER_CMBF_THRESH = 1.8

    # now, each dictionary key is an array (except the month and ifu_amp)
    d = {}

    if not op.exists(cmbf_file):
        d['log'] = "file does not exist"
        return d

    try:
        toks = cmbf_file.split("/")
        month = toks[-2]
        ifu_amp = toks[-1][:7]

        d['month'] = month
        d['ifu_amp'] = ifu_amp

        mflist, mfct = ifu_date_to_multiframe(map_tab, month, toks[-1][1:4])
        # print(f"DEBUG: {month} {toks[-1][1:4]}")

        if mfct > 0:
            # todo: for now just use the first one, if there is more ... or should we use the last??
            # or should we use the one with the longest time for the month
            d['multiframe'] = f"multi_{mflist[0]['ifu_str']}_{ifu_amp[-2:]}"
        else:
            d['multiframe'] = None

        with fits.open(cmbf_file) as hdu:

            if hdu[0].data is not None:
                d['fibnum'] = np.arange(1, len(hdu[0].data)+1, 1)

                # d['sz'] = len(fiber)
                d['mean'] = np.nanmean(hdu[0].data, axis=1)
                d['median'] = np.nanmedian(hdu[0].data, axis=1)
                d['stddev'] = np.nanstd(hdu[0].data, axis=1)

                d['ct_0'] = np.count_nonzero(hdu[0].data == 0, axis=1)
                d['ct_bad'] = np.count_nonzero(hdu[0].data <= BAD_FIBER_CMBF_THRESH, axis=1)
                d['ct_iffy'] = np.count_nonzero(hdu[0].data <= IFFY_FIBER_CMBF_THRESH, axis=1) - d['ct_bad']
                d['ct_high'] = np.count_nonzero(hdu[0].data >= HIGH_FIBER_CMBF_THRESH, axis=1)
                d['status'] = np.full(len(d['fibnum']), -999)  # unset
                d['reason'] = np.full(len(d['fibnum']), 0)  # unset

                d['log'] = "ok"
            else:
                d['fibnum'] = None

                # d['sz'] = len(fiber)
                d['mean'] = None
                d['median'] = None
                d['stddev'] = None

                d['ct_0'] = None
                d['ct_bad'] = None
                d['ct_iffy'] = None
                d['ct_high'] = None
                d['status'] = -999
                d['reason'] = 0
                d['log'] = f"Fail. Empty cmbf file: {cmbf_file}"
    except:
        d['log'] = traceback.format_exc()
        print(f"Exception in get_cmbf_fiber_avg(). {cmbf_file}")

    return d


def evaluate_cmbf_dict(stats_dict):
    """
    what fibers are "bad" or "maybes"?

    fiber 1 (index 0) is at the bottom of DS9 representation

    return -1 = bad, 0 = good, 1 = maybe bad

    """

    if stats_dict is None:
        return stats_dict

    try:
        if len(stats_dict['fibnum']) < 2:
            print(f"Invalid dictionary for evaluate_cmbf_dict. {stats_dict}")
            return stats_dict
    except:
        print(f"Invalid dictionary for evaluate_cmbf_dict. {stats_dict}\n{traceback.format_exc()}")
        return stats_dict

    IFFY_FIBER_CMBF_THRESH = 0.4
    BAD_FIBER_CMBF_THRESH = 0.2
    HIGH_FIBER_CMBF_THRESH = 1.8
    BAD_STDDEV = 1.0
    IFFY_STDDEV = 0.5

    # reasons
    BAD_LOW_MEDIAN = 0x00000001
    BAD_LOW_MEDIAN_CT = 0x00000002
    BAD_STDDEV = 0x00000004
    IFFY_LOW = 0x00000008

    d = stats_dict

    # todo: what logic to establish bad or maybe
    # any with median < BAD is BAD,

    status = np.zeros(len(d['fibnum'])).astype(int)  # reset to 0 = good
    reason = np.zeros(len(d['fibnum'])).astype(int)  # reset to 0 = no reason(s)

    try:
        # first iffy (so it CAN be overwritten by logical order)
        sel = np.array(d['median'] > BAD_FIBER_CMBF_THRESH) * np.array(d['median'] <= IFFY_FIBER_CMBF_THRESH)
        status[sel] = 1
        reason[sel] += IFFY_LOW

        # bad based on median
        sel = np.array(d['median'] <= BAD_FIBER_CMBF_THRESH)
        status[sel] = -1
        reason[sel] += BAD_LOW_MEDIAN

        # bad based on stddev
        sel = np.array(d['stddev'] > BAD_STDDEV)
        status[sel] = -1
        reason[sel] += BAD_LOW_MEDIAN_CT

        # bad based on number of pixels < BAD_FIBER_CMBF_THRESH
        sel = np.array(d['ct_bad'] > 200)  # say 20%?

        d['status'] = status
        d['reason'] = reason
    except:
        print(f"Failure in evaluate_cmbf_dict. dict = {d}\n{traceback.format_exc()}")

    return d


def report_bad_fibers(stats_dict):
    """
    """

    sel = stats_dict['status'] == -1
    print(stats_dict['month'], stats_dict['ifu_amp'], "bad")
    for n in stats_dict['fibnum'][sel]:
        print(f"  {str(n).zfill(3)} {stats_dict['status'][n - 1]} 0x{stats_dict['reason'][n - 1]:08x}")
    print("")


def report_iffy_fibers(stats_dict):
    """
    """

    sel = stats_dict['status'] == 1
    print(stats_dict['month'], stats_dict['ifu_amp'], "iffy")
    for n in stats_dict['fibnum'][sel]:
        print(f"  {str(n).zfill(3)} {stats_dict['status'][n - 1]} 0x{stats_dict['reason'][n - 1]:08x}")
    print("")


def bad_fiber_over_time(stats_dict_array):
    """
    todo: what fibers are constantly bad over time
    """
    pass


#basic (main)
dates = [201701,201702,201703,201704,201705,201706,201707,201708,201709,201710,201711,201712,
         201801,201802,201803,201804,201805,201806,201807,201808,201809,201810,201811,201812,
         201901,201902,201903,201904,201905,201906,201907,201908,201909,201910,201911,201912,
         202001,202002,202003,202004,202005,202006,202007,202008,202009,202010,202011,202012,
         202101,202102,202103,202104,202105,202106,202107,202108,202109,202110,202111,202112,
         202201,202202,202203,202204,202205,202206,202207,202208,202209,202210,202211,202212,
         202301,202302,202303,202304,202305,202306,202307,202308,202309,202310,202311,202312,
         202401,202402,202403,202404,202405,202406,202407]#,202408,202409,202410,202411,202412,
         #202501,202502,202503,202504,202505,202506,202507,202508,202509,202510,202511,202512,
         #202601,202602,202603,202604,202605,202606,202607,202608,202609,202610,202611,202612]

ifuslot = "???"
amp = "??"

all_fT = new_fT()

ifu_map = get_ifu_mapping_table()

cmbfpath = "/scratch/projects/hetdex/lib_calib"
#cmbfpath = "/corral/utexas/Hobby-Eberly-Telesco/lib_calib/"
for date in tqdm(dates):
    date_fT = new_fT()
    fns = sorted(glob.glob(f"{cmbfpath}/{date}/i{str(ifuslot).zfill(3)}a{amp}cmbf.fits"))
    for fn in tqdm(fns):
        d = get_cmbf_all_fiber_stats(fn,ifu_map)
        d = evaluate_cmbf_dict(d)
        fT = new_fT(d)
        date_fT = vstack([date_fT,fT])

    all_fT = vstack([all_fT,date_fT])

all_fT.write("cmbf_stats.dat",format="ascii",overwrite=True)