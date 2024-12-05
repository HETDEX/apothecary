"""

build up the res, sres, and rres*.fits files using the shot.h5 files

USES the full multiframe name, so a spectrograph moved to a new IFU SLOT (for example) will be treated independently

this supersedes the lib_calib/reschi/?res<specid><amp>.fits files
and builds out subdirectories by year+month

by default the calendar time period for the calculation is 5 months as 2 months before and 2 months after the target month
  this can be changed by updating the  prev_months and next_months values

the output is under the CURRENT directory, so you will need to copy to the hetdex directory when done

"""

import numpy as np
import os
import glob
import tables
from astropy.table import Table
import astropy.stats.biweight as biweight
from astropy.io import fits
#from tqdm.notebook import tqdm
from filelock import FileLock
import traceback


basedir_shots = "/scratch/projects/hetdex/hdr5/reduction/data/"
basedir_out = "./lib_calib/reschi/"




#########################
# setup
########################

if not os.path.isdir(basedir_out):
    os.makedirs(basedir_out)
    print(f"Created:  {basedir_out}")

shotfiles_all = sorted(glob.glob(f"{basedir_shots}20*.h5"))
datevshot_all = [os.path.basename(shotfile)[:-3] for shotfile in shotfiles_all]
print(f"Found: {len(shotfiles_all)}")

#read in the list of multiframes and corresponding shots
T_mfs = Table.read("multiframes_unique_slurm.txt",format="ascii")#,names=('multiframe', 'start_shot', 'stop_shot', 'num_shots'))
#update T_mfs with new column for the IFU as multiframe - AMP
#T_mfs['mf_ifu'] = np.array([m[:-3] for m in T_mfs['multiframe']])
T_mfs.sort(['start_shot', 'stop_shot'])
T_mfs['done'] = False

T_badamps = Table.read("/scratch/projects/hetdex/hdr5/survey/amp_flag.fits",format="fits")




#########################
# helpers
#########################

# data table for fiber info
def new_data_table():
    return Table(dtype=[('yyyymm', int),
                        ('shotid', int),
                        ('mf_ifu', (str, 17)),
                        ('multiframe', (str, 20)),
                        ('expnum', int),
                        ('good', int),  # -1 unset, 0 = no, 1 = yes (note backward from bad_amp)
                        ('gid', int)])  # ID to the data
    # ('error1d',(float,(112,1032))),
    # ('sky_subtracted',(float,(112,1032))),
    # ('sky_spectrum',(float,(112,1032))),
    # ])


class FiberData():
    def __init__(self, gid, error1d, sky_subtracted, sky_spectrum):
        # self.yyyymm = yyyymm
        # self.shotid = shotid
        # self.multiframe = multiframe
        # self.expnum = exp_num
        # self.good = good
        self.gid = gid
        self.error1d = error1d  # np.full((112,1032),np.nan)
        self.sky_subtracted = sky_subtracted  # np.full((112,1032),np.nan)
        self.sky_spectrum = sky_spectrum  # np.full((112,1032),np.nan)


class FiberDataDict():

    def __init__(self):
        self.fddict = {}
        self.gid = 0

    def next_gid(self):
        self.gid += 1
        return self.gid

    def reset_gid(self):
        self.gid = 0
        return self.gid

    def insert(self, error1d, sky_subtracted, sky_spectrum):
        gid = self.next_gid()
        fiberdata = FiberData(gid, error1d, sky_subtracted, sky_spectrum)
        self.fddict[gid] = fiberdata
        return gid

    def remove(self, gids):
        try:
            for gid in gids:
                try:
                    del self.fddict[gid]
                except Exception as e:
                    print(f"Exeption in FiberDataDict remove() for gid {gid}: {e}")
        except Exception as e:
            try:
                del self.fddict[gids]  # might be just one and not a list
            except:
                pass
            print(f"Exeption in FiberDataDict remove(): {e}")

    def get_error1d(self, gids):
        try:
            M = []
            for gid in gids:
                try:
                    M.append(self.fddict[gid].error1d)
                except Exception as E:
                    print(f"Could not find {gid} in get_error1d")

            return np.array(M)

        except Exception as e:
            print(f"Exeption in FiberDataDict get_error1d(): {e}")

    def get_sky_subtracted(self, gids):
        try:
            M = []
            for gid in gids:
                try:
                    M.append(self.fddict[gid].sky_subtracted)
                except Exception as E:
                    print(f"Could not find {gid} in get_sky_subtracted")

            return np.array(M)

        except Exception as e:
            print(f"Exeption in FiberDataDict get_sky_subtracted(): {e}")

    def get_sky_spectrum(self, gids):
        try:
            M = []
            for gid in gids:
                try:
                    M.append(self.fddict[gid].sky_spectrum)
                except Exception as E:
                    print(f"Could not find {gid} in get_sky_spectrum")

            return np.array(M)

        except Exception as e:
            print(f"Exeption in FiberDataDict get_sky_spectrum(): {e}")


def is_bad_amp(shotid, multiframe):
    """
    return True if this is a bad_amp (if the flag == 0)
    """
    sel = np.array(T_badamps['shotid'] == shotid) & np.array(T_badamps['multiframe'] == multiframe)
    ct = np.count_nonzero(sel)
    if ct != 1:
        if ct == 0:
            # print(f"bad selection for bad_amps: {shotid} {multiframe} {np.count_nonzero(sel)}")
            return True  # assume bad
        elif ct > 1:
            print(f"bad selection for bad_amps: {shotid} {multiframe} {np.count_nonzero(sel)}")
            return True  # assume bad

    flag = int(T_badamps['flag'][sel])
    if flag:
        return False
    else:
        return True


def get_month_count(t, yyyymm, multiframe):
    try:
        sel = np.array(t['yyyymm'] == yyyymm) & np.array(t['multiframe'] == multiframe)
        ct_all = np.count_nonzero(sel)

        sel = sel & np.array(t['good'] == 1)
        ct_good = np.count_nonzer(sel)

        return ct_good, ct_all

    except Exception as e:
        print(e)
        return -1, -1


# ALL AMPs for an IFU
def insert_data_ifu(t, d, shotid, ifu, h5_fn, skip_bad=True):
    """
    open the shot h5 file for the shotid, read in the 2D arrays and insert into the table
    """
    q_multiframe = None
    try:
        sel = np.array(t['shotid'] == shotid) & np.array(t['mf_ifu'] == ifu)
        if np.count_nonzero(sel) > 0:
            print(f"entry already exists! {shotid} {ifu}")
            return

        # note: this specid, slotid, and ifuid uniquey identify this IFU for THIS shot
        #specid = ifu[6:9]
        q_multiframes = np.array([f'{ifu}_LU', f'{ifu}_LL', f'{ifu}_RU', f'{ifu}_RL'])
        good = np.array([1, 1, 1, 1])  # assume good

        if skip_bad:
            # start_time = time.time()
            for i in range(len(good)):
                if is_bad_amp(shotid, q_multiframes[i]):
                    good[i] = 0
            # print(f"is_bad_amp time: {time.time()-start_time} seconds")
            if np.sum(good) == 0:  # all bad and we are skipping bad, so just bail
                # print(f"skipping bad amp {shotid} {multiframe}")
                return

        query = ""
        for i in range(len(q_multiframes)):
            if skip_bad and good[i] == 0:
                continue
            query += f"(multiframe==b'{q_multiframes[i]}') | "

        if len(query) == 0:
            return
        else:  # strip off the last 3 characters, since this always ends in " | "
            query = query[:-3]

        try:
            # start_time = time.time()
            h5 = tables.open_file(h5_fn)

            exps = h5.root.Data.Fibers.read_where(query, field="expnum")
            mfs = h5.root.Data.Fibers.read_where(query, field="multiframe")

            # wavelengths_all = h5.root.Data.Fibers.read(field="wavelength")
            sky_subtracted = h5.root.Data.Fibers.read_where(query,
                                                            field="sky_subtracted")  # ?? is this the one to use or spectrum or sky_subtracted
            sky_spectrum = h5.root.Data.Fibers.read_where(query,
                                                          field="sky_spectrum")  # ?? is this the one to use or spectrum or sky_subtracted
            # sky_subtracted is "spectrum" - "sky_spectrum" (note this is not calibrated and is in counts)
            error1d = h5.root.Data.Fibers.read_where(query,
                                                     field="error1D")  # rms_all =  h5.root.Data.Fibers.read(field="rms")
            # fiber_ids_all = h5.root.Data.Fibers.read(field="fiber_id") #?? not sure we need this

            sky_subtracted = sky_subtracted.astype(np.float32)
            sky_spectrum = sky_spectrum.astype(np.float32)
            error1d = error1d.astype(np.float32)
            h5.close()
            # print(f"h5 read time: {time.time()-start_time} seconds")

            # start_time = time.time()
            for q_multiframe in q_multiframes:
                for q_exp in [1, 2, 3]:
                    try:
                        sel = np.array(exps == q_exp) & np.array(mfs == q_multiframe.encode())
                        if np.count_nonzero(sel) > 0:
                            # t.add_row([int(str(shotid)[0:6]),shotid,ifu,q_multiframe,q_exp,good[i],error1d[sel],sky_subtracted[sel],sky_spectrum[sel]])
                            gid = d.insert(error1d[sel], sky_subtracted[sel], sky_spectrum[sel])
                            t.add_row([int(str(shotid)[0:6]), shotid, ifu, q_multiframe, q_exp, 1,
                                       gid])  # ,error1d[sel],sky_subtracted[sel],sky_spectrum[sel]])
                    except Exception as e:
                        print('Ex1a', np.shape(error1d[sel]), np.shape(sky_subtracted[sel]),
                              np.shape(sky_spectrum[sel]), e)
            # print(f"tb ({len(t)}) insert time: {time.time()-start_time} seconds")
        except Exception as e:
            print('Ex1', query, shotid, q_multiframe, e)


    except Exception as e:
        print('Ex2', shotid, ifu, e)


def month_add(yyyymm, num):
    year = int(str(yyyymm)[0:4])
    month = int(str(yyyymm)[-2:])

    month = month + num
    if month <= 0:
        while month <= 0:
            month += 12
            year -= 1
    else:
        while month > 12:
            month -= 12
            year += 1

    return int(str(year) + str(month).zfill(2))


def ready_to_compute(t, yyyymm, ifu, prev_months=2, next_months=2):
    """
    return True if we have the data we want to compute res, sres, rres
    (modify as needed)
    assume collection is in date order, so we are not checking the previous months
    """

    try:
        sel = np.array(t['yyyymm'] == month_add(yyyymm, next_months + 1)) & np.array(t['mf_ifu'] == ifu)
        if np.count_nonzero(sel) > 0:
            return True
        else:
            return False
    except Exception as e:
        print(e)
        return False


def reset():
    try:
        del T_data
    except:
        pass


def select_data(t, d, multiframe, month, prev_months=2, next_months=2):
    sel = np.array(t['yyyymm'] >= month_add(month, -1 * prev_months)) & np.array(
        t['yyyymm'] <= month_add(month, next_months))
    sel = sel & np.array(t['multiframe'] == multiframe)

    gids = t['gid'][sel]
    sky_subtracted = d.get_sky_subtracted(gids)
    sky_spectrum = d.get_sky_spectrum(gids)
    error1d = d.get_error1d(gids)

    return gids, sky_subtracted, sky_spectrum, error1d


def compute_all(sky_subtracted, sky_spectrum, error1d):
    bad_value = 0.0
    # todo: need to confirm order of operation ... Matrix Ops first, then biweight?
    error1d[error1d == 0] = np.nan
    sky_spectrum[sky_spectrum == 0] = np.nan
    sky_subtracted[sky_subtracted == 0] = np.nan
    M = sky_subtracted / sky_spectrum
    res = biweight.biweight_location(M, axis=0, ignore_nan=True)
    sres = biweight.biweight_scale(M, axis=0, ignore_nan=True)

    # todo: need to confirm order of operation ... Matrix Ops first, then biweight?
    # there is a very small difference, probabably just precision and is irrelevant since down the 1e-16 range
    # between biweight(sres * M)  vs sres * biweight(M)
    # so, using the latter since it is faster
    M = sky_spectrum / error1d
    rres = sres * biweight.biweight_location(M, axis=0, ignore_nan=True)

    # todo: confirm what to do with nan and zero values?
    # setting to 1.0 seems reasonable based on the consumer (which checks for deviations from 1.0)
    rres[np.isnan(rres)] = bad_value
    rres[rres == 0.0] = bad_value

    sres[np.isnan(sres)] = bad_value
    sres[sres == 0.0] = bad_value

    res[np.isnan(res)] = bad_value
    res[res == 0.0] = bad_value
    return res, sres, rres


def write_status(done, fail, current_mf, current_month):
    try:
        with open("rres_by_ifu_status_slurm.txt", "w") as f:

            f.write(f"current: {current_mf} {current_month}\n\n")

            f.write("Done:\n")
            for k in done.keys():
                f.write(f"{k} {done[k]}\n")

            f.write("\n\nFail:\n")
            for k in fail.keys():
                f.write(f"{k} {fail[k]}\n")

    except Exception as e:
        print(f"write_status failed: {e}")


def write_fits_files(yyyymm, multiframe, res, sres, rres):
    """
    """
    try:
        # make sure path exists
        outdir = os.path.join(basedir_out, str(yyyymm))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
            print(f"Created:  {outdir}")

        # build and write fits to outdirs
        mf = multiframe[6:]
        #amp = multiframe[-2:]

        rres = rres.astype(np.float32)
        sres = sres.astype(np.float32)
        res = res.astype(np.float32)

        # res
        hdu = fits.PrimaryHDU(res)  # essentially with an empty header
        fout = os.path.join(outdir, f'res_{mf}.fits')
        hdu.writeto(fout, overwrite=True)

        # sres
        hdu = fits.PrimaryHDU(sres)  # essentially with an empty header
        fout = os.path.join(outdir, f'sres_{mf}.fits')
        hdu.writeto(fout, overwrite=True)

        # rres
        hdu = fits.PrimaryHDU(rres)  # essentially with an empty header
        fout = os.path.join(outdir, f'rres_{mf}.fits')
        hdu.writeto(fout, overwrite=True)

        return True
    except Exception as e:
        return False

def already_done(ifu):
    done = np.loadtxt("ifu_done.txt",dtype=str)
    if ifu in done:
        return True
    else:
        return False


def get_next_ifu_to_run():
    """
    only run ONE at a TIME ... so multiple instances can only pick up the next IFU as a singelton
    """
    lock = FileLock("rres.lock")
    fn = None
    with lock:
        try:
            next_ifu = sorted(glob.glob("*.open"))[0]
            run_ifu = next_ifu.split(".")[0] + ".run"
            os.rename(next_ifu, run_ifu)
            fn = run_ifu
        except Exception as E:
            print("failed to get next ifu ", E)
            fn = None

    return fn


def mark_ifu_done(fn):
    try:
        done_ifu = fn.split(".")[0] + ".done"
        os.rename(fn, done_ifu)
    except Exception as E:
        print("failed to rename ifu as done ", E)


def parse_ifu_file(fn):
    """
    return the ifu name (without the amp(s)) and the start and stop shots
    """
    try:
        mf_amps, starts, stops, cts, mf_ifus = np.loadtxt(fn, dtype="str", skiprows=1, unpack=True)

        start_shot = np.min(np.array([int(s.replace('v', '')) for s in starts]))
        stop_shot = np.max(np.array([int(s.replace('v', '')) for s in stops]))

        start_shot = str(start_shot)[0:8] + 'v' + str(start_shot)[-3:]
        stop_shot = str(stop_shot)[0:8] + 'v' + str(stop_shot)[-3:]

        return mf_ifus[0], start_shot, stop_shot, mf_amps

    except Exception as E:
        print("failed to parse ifu file ", E)


def find_completed_fits(ifu):
    """
    return list of fits (by month) already complete for the given ifu
    """

    fns = []
    try:
        mf_no_amp = ifu[6:]
        fns = sorted(glob.glob(f"{basedir_out}/*/rres_{mf_no_amp}*.fits"))
        # months = np.unique([f.split("/")[-2] for f in fns])

    except Exception as E:
        pass  # rint("find_completed_fits() ",E)

    return fns


def get_completed_months(completed_list):
    """
    return list of months that are complete or partially complete
    """

    months = None
    try:
        months = np.unique([f.split("/")[-2] for f in completed_list])
    except Exception as E:
        pass

    return list(months)


def get_resume_month(completed_list):
    """
    return target month to resume
    gets list of all completed months and starts at the last month
    which is either complete or partially complete, so worst case, repeats one month
    """

    month = None
    try:
        month = np.unique([f.split("/")[-2] for f in completed_list])[-1]
    except Exception as E:
        pass

    return month


def reset_run():
    """
    mark any *.run files back to *.open
    """
    try:
        run_files = glob.glob("multi*.run")
        for fn in run_files:
            open_fn = fn.split(".")[0] + ".open"
            os.rename(fn, open_fn)
    except Exception as E:
        print("failed reset_run ", E)


def reset_done():
    """
    mark any *.done files back to *.open
    """
    try:
        run_files = glob.glob("multi*.done")
        for fn in run_files:
            open_fn = fn.split(".")[0] + ".open"
            os.rename(fn, open_fn)
    except Exception as E:
        print("failed reset_run ", E)


def check_stop():
    """
    check if the stop file exists ... if so, stop/exit
    :return:
    """
    if os.path.exists("stop"):
        return True
    else:
        return False

##############
# Main loop
#############

# build off of the current month and the calendar adjacent (x) previous and (y) next months
# NOTICE: there might not be the full number of next or previous months as there could be no observations
# or the multiframe is not in use
prev_months = 2
next_months = 2

T_data = new_data_table()
D_data = FiberDataDict()
rres_done = {}
rres_fail = {}

# todo: should we check the output to see if a month+multiframe has already been done?
# caveat: if so, may still need to load the shot.h5 data for that and prior months as they may be used in subsequent months
# so this would have to be a complex check ... see what is the last month completed and the back up prev_months to start
# assuming that there are no in-between months that are missing

# iterating over all multiframes


# todo: change to iterate over all IFUs ... split out the Amps per IFU
mf_ifus = np.array(np.unique(T_mfs['mf_ifu']))  # in roughly date order
ifu_step_size = 10  # use 10 ifus as a time?

all_done = False
mark_done = True #okay to mark done (e.g. this is not a "stop" condition)
next_ifu = None
target_month = None
q_multiframe = None
while not all_done:

    if next_ifu is not None and mark_done:
        mark_ifu_done(next_ifu)

    next_ifu = get_next_ifu_to_run()
    if next_ifu is None:
        all_done = True
        break
    ifu = None
    try:

        ifu, start_shot, stop_shot, q_multiframes = parse_ifu_file(next_ifu)
        months_done = get_completed_months(find_completed_fits(ifu))
        if months_done is None or len(months_done) <= 1:  # new or we resume from the start anyway
            months_done = []  # while techincally could be per amp, just keep as per ifu
            start_idx = datevshot_all.index(start_shot)
            stop_idx = datevshot_all.index(stop_shot)

            for q_multiframe in q_multiframes:
                rres_done[q_multiframe] = []  # initialize a list of months that have been completed (written to rres*fits)
                rres_fail[q_multiframe] = []  # in case there was a problem

            target_month = int(datevshot_all[start_idx].replace('v', '')[0:6])

            print(
                f"Starting {ifu}_XX: {datevshot_all[start_idx]} to {datevshot_all[stop_idx]}, max of {stop_idx - start_idx + 1} shots ...")
        else:  # resume
            last_month = months_done[-1]  # cut off the last one and resume from it
            months_done = months_done[:-1]

            resume_month = month_add(last_month, -1 * prev_months)
            str_resmonth = str(resume_month)
            sd = np.array([d[0:6] == str_resmonth for d in datevshot_all])
            start_shot = np.array(datevshot_all)[sd][0]

            start_idx = datevshot_all.index(start_shot)
            stop_idx = datevshot_all.index(stop_shot)

            print(
                f"Resuming {ifu}_XX: {datevshot_all[start_idx]} to {datevshot_all[stop_idx]}, max of {stop_idx - start_idx + 1} shots ...")

        for datevshot in datevshot_all[start_idx:stop_idx + 1]:

            if check_stop():
                all_done = True
                mark_done = False
                print("Stop file detected. Exiting.")
                break

            shotid = int(datevshot.replace('v', ''))
            if target_month is None:
                target_month = int(datevshot.replace('v', '')[0:6])
                print(f"Reassigning target month {target_month} for {q_multiframe}.")

            # open the shot h5 and read what we want (all 4 amps)
            insert_data_ifu(T_data, D_data, shotid, ifu, os.path.join(basedir_shots, f"{datevshot}.h5"), skip_bad=True)

            # do we have enough to compute a month??
            # if we want 5 months (current +/- 2) what does that look like?

            # for an IFU by month ... all amps
            if ready_to_compute(T_data, target_month, ifu, prev_months=prev_months, next_months=next_months):
                mark_month_done = False
                for q_multiframe in q_multiframes:  # since the amps may not begin/end on the same date, still need to do this by IFU+Amp

                    # do the compute
                    gids, sky_sub, sky_spec, err1d = select_data(T_data, D_data, q_multiframe, target_month)
                    if not ((np.count_nonzero(sky_sub) > 0) and (np.count_nonzero(sky_spec) > 0) and (
                            np.count_nonzero(err1d) > 0)):
                        print(f"No selection to compute for {target_month} {q_multiframe}: ", np.count_nonzero(sky_sub),
                              np.count_nonzero(sky_spec), np.count_nonzero(err1d))
                        continue

                    res, sres, rres = compute_all(sky_sub, sky_spec, err1d)

                    md = np.nanmedian(rres)
                    if md < 0.5 or md > 1.5:
                        print(f"Extreme rres? {md} {q_multiframe} {target_month}")

                    # write out the files
                    files_okay = write_fits_files(target_month, q_multiframe, res, sres, rres)

                    # if the file write was okay:
                    if files_okay:
                        print(f"Wrote {target_month} {q_multiframe} (md: {md:0.2f})")
                        if q_multiframe not in rres_done.keys():
                            rres_done[q_multiframe] = []
                        rres_done[q_multiframe].append(target_month)

                    else:
                        print(f"Failed {target_month} {q_multiframe}")
                        if q_multiframe not in rres_fail.keys():
                            rres_fail[q_multiframe] = []
                        rres_fail[q_multiframe].append(target_month)


                    mark_month_done = True  # at least one completed

                    #write_status(rres_done, rres_fail, q_multiframe, target_month)

                    # remove the old data
                    purge_month = month_add(target_month, -1 * prev_months)
                    sel = np.array(T_data['yyyymm'] == purge_month) & np.array(T_data['multiframe'] == q_multiframe)
                    if np.count_nonzero(sel) > 0:
                        gids = T_data['gid'][sel]
                        T_data.remove_rows(sel)
                        D_data.remove(gids)

                        print(
                            f"Cleaned {np.count_nonzero(sel)} rows for {purge_month} {q_multiframe}. New table size = {len(T_data)}")

                if mark_month_done:
                    months_done.append(target_month)
                    # don't just add 1 to the target_month, there might be gaps, so go to the next possible month, but don't repeat
                    next_target_months = sorted(np.unique(T_data['yyyymm'][T_data['mf_ifu'] == ifu]))
                    idx = 99999
                    if len(next_target_months) > 0:
                        try:
                            if target_month in next_target_months:
                                idx = next_target_months.index(target_month)
                                target_month = next_target_months[idx + 1]
                            else:  # maybe we just deleted (purged) it? could be if there are no more entries for this multiframe
                                # this is the target month we JUST finished
                                max_months = next_target_months[-1]
                                while target_month <= max_months and target_month not in next_target_months:
                                    target_month = month_add(target_month, 1)

                                if target_month > max_months:
                                    target_month = None  # will reset at the top

                        except:
                            target_month = month_add(target_month, 1)


                        while (target_month in months_done) and (idx < len(next_target_months)):
                            idx += 1
                            target_month = next_target_months[idx]

                        if target_month in months_done:
                            # failed to find a new month
                            print(f"Failed to find a new month to run {target_month} {ifu}.")
                            target_month = month_add(target_month, 1)
                    else:  # there are not any more to do
                        print(f"Out of months for {ifu}. Last run = {target_month}")
                        target_month = None

        # end loop over shots for this IFU+amp

        # now do the final stuff
        # need to clear the target_months that did not get run
        for q_multiframe in q_multiframes:
            if check_stop():
                all_done = True
                mark_done = False
                break

            sel = np.array(T_data['multiframe'] == q_multiframe)
            months = np.unique(T_data['yyyymm'][sel])
            print(f"Computing for remaining months ({len(months)})")
            for target_month in months:
                if check_stop():
                    all_done = True
                    mark_done = False
                    break

                gids, sky_sub, sky_spec, err1d = select_data(T_data, D_data, q_multiframe, target_month)
                res, sres, rres = compute_all(sky_sub, sky_spec, err1d)

                md = np.nanmedian(rres)
                if md < 0.5 or md > 1.5:
                    print(f"Extreme rres? {md} {q_multiframe} {target_month}")

                # write out the files
                files_okay = write_fits_files(target_month, q_multiframe, res, sres, rres)

                # if the file write was okay:
                if files_okay:
                    print(f"Wrote {target_month} {q_multiframe} (md: {md:0.2f})")
                    if q_multiframe not in rres_done.keys():
                        rres_done[q_multiframe] = []
                    rres_done[q_multiframe].append(target_month)
                else:
                    if q_multiframe not in rres_fail.keys():
                        rres_fail[q_multiframe] = []
                    rres_fail[q_multiframe].append(target_month)

                #write_status(rres_done, rres_fail, q_multiframe, target_month)

            print(f"Done with {q_multiframe}")

            # remove the old data (don't need by month anymore) we are done with the multiframe, so remove all of it
            sel = np.array(T_data['multiframe'] == q_multiframe)
            if np.count_nonzero(sel) > 0:
                gids = T_data['gid'][sel]
                #T_data = T_data[~sel]
                T_data.remove_rows(sel)
                D_data.remove(gids)


            if q_multiframe in rres_fail.keys() and len(rres_fail[q_multiframe]) > 0:
                print(f"Failures: {rres_fail[q_multiframe]}")

        if mark_done:
            mark_ifu_done(next_ifu)
        next_ifu = None


    except Exception as e:
        print(f"[{ifu}] outer Exception!, {e}")
        try:
            excstr = f"Exception: {e}\n\n{traceback.format_exc()}"
            print(excstr)
            with open(next_ifu+".except","w+") as f:
                f.write(excstr)
                f.write("\n")
        except Exception as e:
            print(f"[{ifu}] outer Exception! Failed to write except file, {e}")



        # try:
        #     rres_fail[row['multiframe']].append(target_month)
        # except:
        #     pass
