"""

used to mark permanently bad fibers in the *cmbf.fits under lib_calib
this is by fiber NUMBER (1-112) for a single IFU(Slot) + Amp

the original file is backedup locally, then the cmbf.fits file is updated with the specified fibers being set to zero

"""
import numpy as np
import os
import sys
import glob
import astropy.io.fits as fits
import shutil
import traceback
import copy
from datetime import datetime

#cmbf_topdir = "/scratch/03261/polonius/bad_fibers/cmbf_test/lib_calib"
cmbf_topdir = "/scratch/projects/hetdex/lib_calib"



args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

NO_ACTION = False #for sompatibility with zero_all_cmbf_for_month

#parm order:   IFU(slot), AMP, comma separated fiber numbers, startdate, enddate, optional path

#the first five are positional and required
# any exceptions are fatal

ifuslot = int(args[0])
amp = str(args[1])

#comma list
toks = args[2].split(",")
fiber_number_list = [int(t) for t in toks]

startmonth = int(args[3])
stopmonth = int(args[4])

if len(args) > 5:
    cmbf_topdir = args[5]



def month_range(startmonth, stopmonth):
    """All months from a to b inclusive, as YYYYMM ints, ascending."""

    def to_ord(ym):
        y, m = divmod(ym, 100)
        if not 1 <= m <= 12:
            raise ValueError(f"invalid month in {ym}")
        return y * 12 + (m - 1)

    lo, hi = sorted((startmonth, stopmonth))
    return [(o // 12) * 100 + (o % 12) + 1 for o in range(to_ord(lo), to_ord(hi) + 1)]


def get_cmbf_file_list(ifu, amp, startmonth=None, stopmonth=None, topdir=None):
    """
    ifu is IFUSlot
    dates are inclusive

    returns a list

    """
    try:

        cmbf_filename = f"i{str(int(ifu)).zfill(3)}a{amp}cmbf.fits"
        search_list = []
        fns = []

        if topdir is None:
            topdir = cmbf_topdir

        if startmonth is not None and stopmonth is not None:
            if int(stopmonth) < int(startmonth):
                print("Invalid date range")
                return []
            # make a list of all months in range
            month_list = month_range(startmonth, stopmonth)
            for month in month_list:
                search_list.append(os.path.join(topdir, str(month), cmbf_filename))
        else:
            search_list = [os.path.join(topdir, cmbf_filename)]

        for search_dir in search_list:
            fns.append(glob.glob(search_dir))

        fns = sum(fns, [])  # fine unless the list size gets really big (which should not happen given the
        # limitation of months of recent years)
    except:
        print(traceback.format_exc())

    return fns


# def zero_fibers_one_file(fn, fiber_number_list=[]):
#     """
#
#     fibers index from the bottom
#
#     """
#     try:
#         # first make a backup
#
#         if len(fiber_number_list) == 0 or len(fiber_number_list) > 111:
#             print("Invalid fiber_number_list. Bad length.")
#             return -1
#         elif np.any(np.array(fiber_number_list) < 1) or np.any(np.array(fiber_number_list) > 112):
#             print("Invalid fiber_number_list. Out of range values.")
#             return -1
#
#         if os.path.exists(fn):
#             fnb = fn.replace("cmbf.fits", "cmbf.backup.fits")
#             target_ready = False
#             ct = 0
#             while not target_ready:
#                 if os.path.exists(fnb):
#                     ct += 1
#                     fnb = fn.replace("cmbf.fits", f"cmbf.backup_{ct}.fits")
#                 else:
#                     target_ready = True
#             shutil.copy(fn, fnb)
#             print(f"Backup: {fnb}")
#         else:
#             print("File does not exist")
#             return -1
#
#         with fits.open(fn, mode='update') as hdu:
#             for fnum in fiber_number_list:
#                 hdu[0].data[fnum - 1] *= 0
#
#         print(f"Updated: {fn}\n")
#
#     except:
#         print(traceback.format_exc())
#         return -1
#
#     return 0

def zero_fibers_one_file(fn, fiber_number_list=[]):
    """

    fibers index from the bottom

    """
    try:
        # first make a backup

        if len(fiber_number_list) == 0 or len(fiber_number_list) > 111:
            print(f"{fn} Invalid fiber_number_list. Bad length.")
            return -1
        elif np.any(np.array(fiber_number_list) < 1) or np.any(np.array(fiber_number_list) > 112):
            print(f"{fn} Invalid fiber_number_list. Out of range values.")
            return -1

        do_change = False
        if os.path.exists(fn):

            #is a change going to happen? could be the current .fits is already marked
            with fits.open(fn, mode='readonly') as hdu:
                data = copy.copy(hdu[0].data)
                for fnum in fiber_number_list:
                    data[fnum - 1] *= 0

                if np.any(data - hdu[0].data):
                    #there is going to be a change
                    do_change = True

                if not NO_ACTION:  # nice double negative
                    if do_change:
                        #notice: the largest number is then the most recent backup
                        fnb = fn.replace("cmbf.fits", "cmbf.backup_1.fits")
                        target_ready = False
                        ct = 0
                        while not target_ready:
                            if os.path.exists(fnb):
                                ct += 1
                                fnb = fn.replace("cmbf.fits", f"cmbf.backup_{ct}.fits")
                            else:
                                target_ready = True
                        shutil.copy(fn, fnb)
                        print(f"Backup: {fnb}")

                        with fits.open(fn, mode='update') as hdu:
                            now = datetime.now().strftime("%Y-%m-%d")
                            for fnum in fiber_number_list:
                                hdu[0].data[fnum - 1] *= 0
                                hdu[0].header.add_history(f"{now} set fiber number {str(fnum).zfill(3)} to zero. Bad fiber.")

                        print(f"{fn} updated. Fibers: {fiber_number_list}")
                    else:
                        #DEFBUG
                        print(f"{fn} already updated.")
                else:
                    if do_change:
                        print(f"{fn} would be updated. Fibers: {fiber_number_list}")

        else:
            print(f"{fn} File does not exist")
            return -1

    except:
        print(traceback.format_exc())
        return -1

def zero_fibers(ifu, amp, fiber_number_list=[], startmonth=None, stopmonth=None, topdir=None):
    """

    """

    try:

        file_list = get_cmbf_file_list(ifu, amp, startmonth, stopmonth, topdir)
        for fn in file_list:
            try:
                zero_fibers_one_file(fn, fiber_number_list)
            except:
                print(f"Failed: {fn}\n", traceback.format_exc())

    except:
        print(traceback.format_exc())
        return -1



def main():

    zero_fibers(ifuslot,amp,fiber_number_list,startmonth,stopmonth,cmbf_topdir)

if __name__ == "__main__":
    main()