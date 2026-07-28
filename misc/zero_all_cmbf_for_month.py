"""

used to mark permanently bad fibers in the *cmbf.fits under lib_calib
this is by fiber NUMBER (1-112) for a single IFU(Slot) + Amp

the original file is backedup locally, then the cmbf.fits file is updated with the specified fibers being set to zero


!!! Based on zero_cmbf.py   BUT this one uses the HETDEX_API badfib.tab file and operates on a single month
    of lib_calib !!!


    This version takes no parameters and must be run within the lib_calib/yyyymm  that is the target


    Needs
    HETDEX_API/hetdex_api/known_issues/hdr5/badfib.tab
    apothecary/misc/ifu_mf_map_full_hetdex.dat (so can map the full multiframe to the IFU slot for the lib_calib month)



"""
import copy

import numpy as np
import os
import sys
import glob
import astropy.io.fits as fits
import shutil
import traceback
from astropy.table import Table
from datetime import datetime

#import hetdex_api
#hetdex_api_path = hetdex_api.__path__

badfib_path = "/corral-repl/utexas/Hobby-Eberly-Telesco/hdrX/software/hetdex_api/known_issues/hdr5/badfib.tab"
ifu_map_path = "/corral-repl/utexas/Hobby-Eberly-Telesco/hdrX/software/apothecary/misc/ifu_mf_map_full_hetdex.dat"

#cmbf_topdir = "/scratch/03261/polonius/bad_fibers/cmbf_test/lib_calib"
cmbf_topdir = "/scratch/projects/hetdex/lib_calib"



args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

NO_ACTION =True
if args[0] == "-n":
    #no action ... list list what would be modified
    NO_ACTION = True
elif args[0] == "-x":
    NO_ACTION = False
else:
    print("Safety ... did not specify -n or -x")
    exit(0)


#parm order:   IFU(slot), AMP, comma separated fiber numbers, startdate, enddate, optional path

#the first five are positional and required
# any exceptions are fatal

# ifuslot = int(args[0])
# amp = str(args[1])
#
# #comma list
# toks = args[2].split(",")
# fiber_number_list = [int(t) for t in toks]
#
# startmonth = int(args[3])
# stopmonth = int(args[4])
#
# if len(args) > 5:
#     cmbf_topdir = args[5]



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
                                if np.any(hdu[0].data[fnum - 1]):
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

    return 0


# def zero_fibers(ifu, amp, fiber_number_list=[], startmonth=None, stopmonth=None, topdir=None):
#     """
#
#     """
#
#     try:
#
#         file_list = get_cmbf_file_list(ifu, amp, startmonth, stopmonth, topdir)
#         for fn in file_list:
#             try:
#                 zero_fibers_one_file(fn, fiber_number_list)
#             except:
#                 print(f"Failed: {fn}\n", traceback.format_exc())
#
#     except:
#         print(traceback.format_exc())
#         return -1
#


def main():

    #assume we are under <path>/lib_calib/yyyymm
    cwd = os.getcwd()
    current_yyyymm = int(os.path.basename(cwd)) #if this is NOT where we are, it will throw an exception and abend
    if not 201701 <= current_yyyymm < 20990101:
        print(f"Fail. Nonsense yyyymm {current_yyyymm}.")
        exit(-1)


    cmbf_fns = sorted(glob.glob("i???a??cmbf.fits"))

    if len(cmbf_fns) == 0:
        print("Fail. No *cmbf.fits files found in the current directory.")
        exit(-1)

    IFUMAP = Table.read(ifu_map_path,format="ascii") #want most values as ints
    IFUMAP_last_yyyymm = sorted(np.unique(IFUMAP['yyyymm']))[-1]
    BADFIB = Table.read(badfib_path,format="ascii")


    if current_yyyymm > IFUMAP_last_yyyymm:
        lookup_yyyymm = IFUMAP_last_yyyymm
    else:
        lookup_yyyymm = current_yyyymm

    for fn in cmbf_fns:
        #does this have an appropriate entry in the IFUMAP table for this date?
        #assume the last mapping for any date AFTER the end of the table
        # yyyymmdd yyyymm specid ifuslot ifuid (all read in as ints)

        bn = os.path.basename(fn)
        ifuslot = int(bn[1:4]) #iXXXaXXcmbf.fits
        specid = 0
        ifuid = 0
        amp = bn[5:7]
        sel = np.array(IFUMAP['yyyymm']==lookup_yyyymm) * np.array(IFUMAP['ifuslot']==ifuslot)
        #should be 1 or 0 hits
        ct = np.count_nonzero(sel)
        if ct == 0:
            print(f"Problem. {ct} hits for {fn}")
            continue
        else: #this is usually okay as there are many days in a month
            if len(np.unique(IFUMAP['specid'][sel])) != 1 or len(np.unique(IFUMAP['ifuid'][sel])) != 1:
                #now we have a problem
                print(f"Problem. {ct} hits for {fn}. Could not uniquely identify the multiframe.")
                print(f"    specid={np.unique(IFUMAP['specid'][sel])} , ifuid={np.unique(IFUMAP['ifuid'][sel])}")
                continue
            else:
                #all good
                specid = IFUMAP['specid'][sel][0]
                ifuid = IFUMAP['ifuid'][sel][0]

        mf = f"multi_{str(specid).zfill(3)}_{str(ifuslot).zfill(3)}_{str(ifuid).zfill(3)}_{amp}"

        #now, find out if there are any hits in the BADFIB
        sel = BADFIB['multiframe']==mf  #there can be 0 or more
        ct = np.count_nonzero(sel)

        if ct > 0:
            #print("")
            zero_fibers_one_file(fn, sorted(list(BADFIB['fibnum'][sel])))

if __name__ == "__main__":
    main()