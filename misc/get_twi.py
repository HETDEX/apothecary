
import sys
import os
import glob
import tarfile
from tqdm import tqdm



args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

#first two parms are unlabled

if len(args) > 0:
    yyyymmdd = int(args[0])

if len(args) > 1:
    ifu_pattern = args[1]

do_extract = False
if len(args) > 2:
    if "extract" in args[2]:
        do_extract = True

prefix = None
if len(args) > 3:
    prefix = args[3]

first_only = False
if len(args) > 4:
    if "first" in args[4]:
        first_only = True

top_path = "/work/03946/hetdex/maverick/" #first choice, use maverick
alt_top_path = "/corral-repl/utexas/Hobby-Eberly-Telesco/het_raw"
use_dir = False
#top_path = "/scratch/03261/polonius/bad_fibers"

def check_basepath(yyyymmdd):
    """

    :param yyyymdd:
    :return:
    """

    global top_path, alt_top_path, use_dir

    if os.path.exists(os.path.join(top_path,f"{yyyymmdd}")):
        #all good, do nothing
        use_dir = True
        return top_path
    elif os.path.exists(os.path.join(alt_top_path,f"{yyyymmdd}.tar")):
        #update and return
        top_path = alt_top_path
        use_dir = False
        return top_path
    elif os.path.exists(os.path.join(alt_top_path,f"{yyyymmdd}.tar")):
        #the directory is here, not the tar
        top_path = alt_top_path
        use_dir = False
        return top_path
    else:
        print("Cannot find any paths")
        return None

def find_nested_twi_from_tar(yyyymmdd, ifu_pattern=None, do_extract=False, first=False):
    """
    todo: maybe allow an IFU or IFU+amp (ifu_pattern could be 096 or 096R or 096RU)
    top level is still a tar, e.g. /corral/utexas/Hobby-Eberly-Telesco/het_raw/

    """

    all_twilights = []
    all_outfiles = []

    if type(yyyymmdd) is int:
        yyyymmdd = str(yyyymmdd)

    # toplevel tar
    fn = os.path.join(top_path, f"{yyyymmdd}.tar")
    if os.path.exists(fn):
        with tarfile.open(fn, 'r') as toptar:
            tar_to_check = []
            topname = toptar.next()
            while topname is not None:
                # print(f"DEBUG: {topname.path}, type = {type(topname)}")
                topname = topname.path
                if "virus0" in topname:
                    tar_to_check.append(topname)
                topname = toptar.next()

            done_done = False
            for topname in tqdm(tar_to_check):
            #for topname in tar_to_check:
                if done_done:
                    break
                # print(f"DEBUG: opening from tar: {topname}")
                twilights = []
                outfiles = []
                with toptar.extractfile(topname) as subtar:
                    with tarfile.open(fileobj=subtar, mode="r") as tarf:
                        sub_done = False
                        while not sub_done:
                            nextname = tarf.next()
                            if nextname is None:
                                sub_done = True
                                # print(f"sub_done")
                                # print(f"toptar(2)={toptar}")
                            else:
                                nextname = nextname.path
                                #print(f"DEBUG: checking: {nextname}")
                                if nextname[-8:] == "twi.fits":
                                    found = False
                                    if ifu_pattern is not None:
                                        if ifu_pattern in os.path.basename(nextname):
                                            found = True
                                            #print(f"Found: {nextname}")
                                    else:
                                        found = True
                                        #print(f"Found: {nextname}")

                                    if found:
                                        twilights.append(nextname)
                                        if first:
                                            sub_done = True
                                            done_done = True
                                            break
                                else: #these might be science frames or dark frames or bias, etc; no point in continuuing
                                    sub_done = True
                                    break

                        #doing the extract must also advance the file pointer, so it messes up
                        #the *.next() call
                        if do_extract:
                            for nextname in twilights:
                                with tarf.extractfile(nextname) as src:
                                    outname = os.path.basename(nextname)
                                    if prefix is not None:
                                        outname = f"{prefix}_{outname}"
                                    with open(outname, "wb") as dest:
                                        dest.write(src.read())
                                        outfiles.append(outname)
                        all_twilights += twilights
                        all_outfiles += outfiles
    else:
        print("path does not exist")

    #print(len(twilights),len(outfiles))
    return all_twilights, all_outfiles


def find_nested_twi_from_dir(yyyymmdd, ifu_pattern=None, do_extract=False, first=False):
    """
    todo: maybe allow an IFU or IFU+amp (ifu_pattern could be 096 or 096R or 096RU)
    top level is a directory, e.g. /work/03946/hetdex/maverick/YYYYMMDD

    """

    all_twilights = []
    all_outfiles = []

    if type(yyyymmdd) is int:
        yyyymmdd = str(yyyymmdd)

    # toplevel tar
    pn = os.path.join(top_path, f"{yyyymmdd}")
    if os.path.exists(pn):
        tar_to_check = sorted(glob.glob(f"{pn}/virus/virus0*.tar"))

        if len(tar_to_check) > 0:
            done_done = False
            for subtar in tqdm(tar_to_check):
                if done_done:
                    break
                # print(f"DEBUG: opening from tar: {topname}")
                twilights = []
                outfiles = []
                if True: #just here to keep the same indentation as the find_nested_twi_from_tar() function
                    with tarfile.open(subtar, mode="r") as tarf:
                        sub_done = False
                        while not sub_done:
                            nextname = tarf.next()
                            if nextname is None:
                                sub_done = True
                                # print(f"sub_done")
                                # print(f"toptar(2)={toptar}")
                            else:
                                nextname = nextname.path
                                #print(f"DEBUG: checking: {nextname}")
                                if nextname[-8:] == "twi.fits":
                                    found = False
                                    if ifu_pattern is not None:
                                        if ifu_pattern in os.path.basename(nextname):
                                            found = True
                                            #print(f"Found: {nextname}")
                                    else:
                                        found = True
                                        #print(f"Found: {nextname}")

                                    if found:
                                        twilights.append(nextname)
                                        if first:
                                            sub_done = True
                                            done_done = True
                                            break
                                else: #these might be science frames or dark frames or bias, etc; no point in continuuing
                                    sub_done = True
                                    break

                        #doing the extract must also advance the file pointer, so it messes up
                        #the *.next() call
                        if do_extract:
                            for nextname in twilights:
                                with tarf.extractfile(nextname) as src:
                                    outname = os.path.basename(nextname)
                                    if prefix is not None:
                                        outname = f"{prefix}_{outname}"
                                    with open(outname, "wb") as dest:
                                        dest.write(src.read())
                                        outfiles.append(outname)
                        all_twilights += twilights
                        all_outfiles += outfiles
    else:
        print("path does not exist")

    #print(len(twilights),len(outfiles))
    return all_twilights, all_outfiles


def main():

    if yyyymmdd < 210000: #probably just a YYYYMM
        to_check_tars = os.path.join(top_path,f"{yyyymmdd}??.tar")
        print(f"Listing dates to check for month: {yyyymmdd} ... ")
        fns1 = sorted(glob.glob(os.path.join(top_path,f"{yyyymmdd}")))
        fns2 = sorted(glob.glob(os.path.join(alt_top_path,f"{yyyymmdd}??.tar")))
        fns = list(fns1) + list(fns2)
        for f in fns:
            print(os.path.basename(f))
    else:
        check_basepath(yyyymmdd)

        print(f"Using {top_path}")
        if use_dir:
            twi, out = find_nested_twi_from_dir(yyyymmdd, ifu_pattern=ifu_pattern, do_extract=do_extract, first=first_only)
        else:
            twi,out = find_nested_twi_from_tar(yyyymmdd, ifu_pattern=ifu_pattern, do_extract=do_extract,first=first_only)
        if len(out) > 0:
            for t,o in zip(twi,out):
                print(t,"-->",o)
        else:
            for t in twi:
                print(t)

if __name__ == "__main__":
    main()