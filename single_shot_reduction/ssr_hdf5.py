"""
Single Shot Reduction HDF5 file creation

For long term storage of key data for a single shot reduction

Based on the elixer_hdf5.py code AND the HETDEX_API shot h5 code

This file is for a single shot (observation) ONLY. Do NOT commbine shots.

#note: for vm-small, 3 shots at a time seems to be about the memory limit

"""

# 0.1.0 long-running baseline without notes
# 0.1.1 added FullSkyModel groups with reduced precision from float64 to float32
# 0.1.2 restored VIRUSImage data (the CCD clean_image, image(with sky) and error, but really try to keep as float16
#       the images can force to float32 but if only the error is large, cap the error to max float16
# 0.1.3 restored FiberIndex as its own separate table to speed up searches
# 0.1.4 add in the Astrometry Collapsed images as they can come in handy for debugging and are part of analysis notebook
#       turned on group level compression for these iamges AND for the CCD images
#       added in Calibration.Throughput.throughput
# 0.1.5 repace --exclude_ccd_images with --minimum as applies to more than just the VIRUS CCD images
#       include all additional original hdf5 info

__version__ = '0.1.5'


import numpy as np
import tables
import os
from pathlib import Path
import sys
import glob
from PIL import Image
from tqdm import tqdm
import traceback
import time

try:
    from filelock import FileLock
except:
    print("You need to install filelock (e.g.: pip install --user filelock) ")
    exit(-1)


UNSET_FLOAT = -999.999
UNSET_INT = -99999
UNSET_STR = ""
UNSET_NAN = np.nan
SUPPORTED_ELIXER_H5_VERSIONS = [b"0.9.2",b"0.10.0"]
SHOW_TQDM = True

Include_CoaddImage_matched = True
Include_CoaddImage_ifugrid = True
Include_CoaddImage_skyview = True

HalfFloatCol = tables.Float16Col
FullFloatCol = tables.Float32Col
DblFloatCol = tables.Float64Col
StdFloatCol = HalfFloatCol #tables.Float16Col #could change to float32 if needed


HalfFloat = np.float16
FullFloat = np.float32
DblFloat = np.float64
StdFloat = HalfFloat #np.float16  #could change to float32 if needed

threshold_16 = 6e4     #if abs value greater than this, go with float32 (just due to precision loss), about half-max
                       #float16 can sort of represent up to about 65K (65504.0 is max)
                       #AND be careful about where to apply this. A loss of precision is fine for fluxes where the
                       #uncertianties are typically around 20-25%, but it is NOT ok for, say, wavelengths which have
                       #uncertainties way less than 1%

image_threshold_16 = 64800.0 #back off the 65504 just a bit. experimentally, I like this value for rounding better
clip_float16 = 65500.0

FATAL_EXIT = 0 #global flag to terminate
elixer_h5_force = False  #if True (set later) ignore the version constraint and try anyway
#exclude_ccd_images = False #do not include the CCD images (VIRUSImage table)
minimum_h5 = False #replaced exclude_ccd_images
full_h5 = True

#sanity force values
if full_h5:
    minimum_h5 = False
    Include_CoaddImage_matched = True
    Include_CoaddImage_ifugrid = True
    Include_CoaddImage_skyview = True

#logging will just be prints
#this is all single threaded, no real management needed
class PrintLog:
    def __init__(self,logfile=None):
        self.logfile =logfile

    def log(self,msg,exc_info=False):
        print(msg,flush=True)
        if self.logfile:
            try:
                with open(self.logfile,"a") as f:
                    if exc_info:
                        f.write(msg + traceback.format_exc() + "\n")
                    else:
                        f.write(msg+"\n")
            except:
                pass

    def critical(self,msg,exc_info=False):
        self.log(msg,exc_info)

    def error(self,msg,exc_info=False):
        self.log(msg,exc_info)

    def info(self,msg,exc_info=False):
        self.log(msg,exc_info)

    def warning(self,msg,exc_info=False):
        self.log(msg,exc_info)

    def debug(self,msg,exc_info=False):
        self.log(msg,exc_info)
# end class def

log = PrintLog('h5_merge.log')

Lock_tmp_mutex_fn = "tmp_ssrcompress.mutex"
Lock_tmp_ct_fn = "tmp_ssrcompress.ct"
Max_Simultaneous_Shots = 3 #this depends on options

SafeActiveShotsSleep = 30.0 # recheck every xx seconds

def wait_to_run(max_procs=3,datevshot="???",clean_up=False): #,safelimit=0):
    """

    use filelock and counts to limit the number of SSR build running simultaneously

    :param max_procs: how many can run at time

    """
    sleep_secs = 1.0 if clean_up else SafeActiveShotsSleep

    try:
        lock = FileLock(Lock_tmp_mutex_fn)
        abort = False

        if max_procs > 0:
            redlight = True
            if clean_up:
                print(f"[{datevshot}] cleaning up ...")
            else:
                print(f"[{datevshot}] checking if safe to start ...")

            while redlight:
                with lock:
                    #how many are active?
                    with open(Lock_tmp_ct_fn,"a+") as f:
                        f.seek(0)
                        ct = f.readline()
                        try:
                            ct = int(ct)
                        except:
                            ct = 0

                        if ct < 0 and not clean_up: #this is the signal to abort and exit
                            print(f"[{datevshot}] recieved mutex count abort. Exiting ...")
                            abort = True


                        f.seek(0)
                        if clean_up: #we are done and need to remove THIS runner
                            ct = max(0,ct-1)
                            f.truncate()
                            f.write(f"{ct}\n")
                            redlight = False
                        else:
                            if ct < max_procs and not abort:
                                ct +=1
                                f.truncate()
                                f.write(f"{ct}\n")
                                redlight = False
                            else:
                                #still need to wait
                                pass


                if not redlight:
                    if clean_up:
                        print(f"[{datevshot}] Done.")
                    elif not abort:
                        print(f"[{datevshot}] cleared to start.")
                    else:
                        exit(-1)
                else:
                    #print(f"[{datevshot}] too many active shots. Must wait ...")
                    time.sleep(sleep_secs)
        # lock auto releases
    except:
        print(f"[{datevshot}] Exception! in wait_to_run()",traceback.format_exc())

#make a class for each table
class Version(tables.IsDescription):
    #version table, very basic info
    version = tables.StringCol(itemsize=16, dflt='',pos=0)
    version_pytables = tables.StringCol(itemsize=16, dflt='',pos=1)


class Detections(tables.IsDescription):
    """
    mostly a clone of elixer's h5 Detections table, with the shot specific info removed as redundant to the shot table
    """

    detectid = tables.Int64Col(pos=0)
    elixer_version = tables.StringCol(itemsize=16,pos=1) #version of elixer that generated this detection report
    elixer_datetime = tables.StringCol(itemsize=21,pos=2) #YYYY-MM-DD hh:mm:ss


    # h5_report_idx = tables.Int32Col(pos=3,dflt=-1)
    # h5_report_gpcd = tables.StringCol(itemsize=32,pos=4,dflt="")
    # h5_neighbor_idx = tables.Int32Col(pos=5,dflt=-1)
    # h5_neighbor_gpcd = tables.StringCol(itemsize=32,pos=6,dflt="")

    h5_report_idx = tables.Int32Col(pos=3,dflt=-1)
    h5_report_id = tables.Int32Col(pos=4,dflt=-1)
    h5_neighbor_idx = tables.Int32Col(pos=5,dflt=-1)
    h5_neighbor_id = tables.Int32Col(pos=6,dflt=-1)

    #shotid = tables.Int64Col() #redundant with shot table
    #obsid = tables.Int32Col() #redundant with shot table

    specid = tables.StringCol(itemsize=3)
    ifuslot = tables.StringCol(itemsize=3)
    ifuid = tables.StringCol(itemsize=3)
    amp = tables.StringCol(itemsize=2)
    multiframe = tables.StringCol(itemsize=20)

    # seeing_fwhm = tables.Float32Col(dflt=UNSET_FLOAT) #redundant with shot table
    # response = tables.Float32Col(dflt=UNSET_FLOAT) #redundant with shot table
    # fieldname = tables.StringCol(itemsize=32) #redundant with shot table


    #about the detection

    ra = tables.Float32Col(dflt=UNSET_FLOAT,pos=7)
    dec = tables.Float32Col(dflt=UNSET_FLOAT,pos=8)
    wavelength_obs = tables.Float32Col(dflt=UNSET_FLOAT,pos=9)
    wavelength_obs_err = tables.Float32Col(dflt=UNSET_FLOAT,pos=10)
    apcor_4500 = tables.Float32Col(dflt=UNSET_FLOAT)

    z_best = tables.Float32Col(dflt=-1.0,pos=11)
    z_best_pz = tables.Float32Col(dflt=0.0,pos=12)


    z_best_plya_thresh = tables.Float32Col(dflt=-1.0, pos=13)
    z_best_2 = tables.Float32Col(dflt=-1.0, pos=14)
    z_best_pz_2 = tables.Float32Col(dflt=0.0, pos=15)
    z_best_plya_thresh_2 = tables.Float32Col(dflt=-1.0, pos=16)
    z_best_3 = tables.Float32Col(dflt=-1.0, pos=17)
    z_best_pz_3 = tables.Float32Col(dflt=0.0, pos=18)
    z_best_plya_thresh_3 = tables.Float32Col(dflt=-1.0, pos=19)


    flags = tables.Int32Col(dflt=0,pos=20)
    review = tables.Int8Col(dflt=0,pos=21)
    cluster_parent = tables.Int64Col(dflt=0,pos=22)


    flux_line = tables.Float32Col(dflt=UNSET_FLOAT) #actual flux not flux density
    flux_line_err = tables.Float32Col(dflt=UNSET_FLOAT)

    #new 0.9.0
    flux_line_obs = tables.Float32Col(dflt=UNSET_FLOAT)  # not dustcorrected
    flux_line_obs_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_line_dust_corr = tables.Float32Col(dflt=UNSET_FLOAT)

    fwhm_line_aa = tables.Float32Col(dflt=UNSET_FLOAT)
    fwhm_line_aa_err = tables.Float32Col(dflt=UNSET_FLOAT)
    fwhm_line_kms = tables.Float32Col(dflt=UNSET_FLOAT)
    fwhm_line_kms_err = tables.Float32Col(dflt=UNSET_FLOAT)
    sn = tables.Float32Col(dflt=UNSET_FLOAT)
    sn_err = tables.Float32Col(dflt=UNSET_FLOAT)
    chi2 = tables.Float32Col(dflt=UNSET_FLOAT)
    chi2_err = tables.Float32Col(dflt=UNSET_FLOAT)

    continuum_line = tables.Float32Col(dflt=UNSET_FLOAT) #continuum (y-offset) from Gaussian fit to the line
    continuum_line_err = tables.Float32Col(dflt=UNSET_FLOAT)

    #new 0.9.0
    continuum_line_obs = tables.Float32Col(dflt=UNSET_FLOAT)  # not dustcorrected
    continuum_line_obs_err = tables.Float32Col(dflt=UNSET_FLOAT)

    eqw_rest_lya_line = tables.Float32Col(dflt=UNSET_FLOAT)
    eqw_rest_lya_line_err = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_line = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_line_max = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_line_min = tables.Float32Col(dflt=UNSET_FLOAT)


    continuum_wide = tables.Float32Col(dflt=UNSET_FLOAT) #continuum from (mostly) the full spectrum width, from pseudo g-band magnitude
    continuum_wide_err = tables.Float32Col(dflt=UNSET_FLOAT)
    mag_g_wide = tables.Float32Col(dflt=UNSET_FLOAT) #the pseudo g-band magnitude from the spectrum
    mag_g_wide_err = tables.Float32Col(dflt=UNSET_FLOAT)
    mag_g_wide_limit = tables.Float32Col(dflt=24.8)
    eqw_rest_lya_wide = tables.Float32Col(dflt=UNSET_FLOAT)
    eqw_rest_lya_wide_err = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_wide = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_wide_max = tables.Float32Col(dflt=UNSET_FLOAT)
    plae_wide_min = tables.Float32Col(dflt=UNSET_FLOAT)

    continuum_masked = tables.Float32Col(dflt=UNSET_FLOAT) #continuum from (mostly) the full spectrum width, from pseudo g-band magnitude
    continuum_masked_err = tables.Float32Col(dflt=UNSET_FLOAT)


    #ELiXer solution based on extra lines
    multiline_flag = tables.BoolCol(dflt=False) #True if s a single "good" solution
    multiline_z = tables.Float32Col(dflt=UNSET_FLOAT)
    multiline_rest_w = tables.Float32Col(dflt=UNSET_FLOAT)
    multiline_prob = tables.Float32Col(dflt=UNSET_FLOAT)
    multiline_raw_score = tables.Float32Col(dflt=UNSET_FLOAT)
    multiline_frac_score = tables.Float32Col(dflt=UNSET_FLOAT)
    multiline_name = tables.StringCol(itemsize=16)

    pseudo_color_flag = tables.Int64Col(dflt=0)

    pseudo_color_blue_flux = tables.Float32Col(dflt=UNSET_FLOAT) #all un uJy
    pseudo_color_blue_flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    pseudo_color_red_flux = tables.Float32Col(dflt=UNSET_FLOAT)
    pseudo_color_red_flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    pseudo_color_rvb_ratio = tables.Float32Col(dflt=UNSET_FLOAT)
    pseudo_color_rvb_ratio_err = tables.Float32Col(dflt=UNSET_FLOAT)

    #ELiXer combined (rules, inv variance, weights and Bayes) classification info
    combined_plae = tables.Float32Col(dflt=UNSET_FLOAT)   #combination of all PLAE/POII
    combined_plae_lo = tables.Float32Col(dflt=UNSET_FLOAT)
    combined_plae_hi = tables.Float32Col(dflt=UNSET_FLOAT)

    plya_classification = tables.Float32Col(dflt=UNSET_FLOAT) #final, combine P(LAE) (0.0 - 1.0) #fomerly plae_classification

    spurious_reason = tables.StringCol(itemsize=32,dflt=UNSET_STR)
    combined_continuum = tables.Float32Col(dflt=UNSET_FLOAT)   #combination of all continuum estimates
    combined_continuum_err = tables.Float32Col(dflt=UNSET_FLOAT)
    combined_continuum_err_stat = tables.Float32Col(dflt=UNSET_FLOAT)
    combined_continuum_nondetect = tables.Int32Col(dflt=UNSET_INT)


    combined_eqw_rest_lya = tables.Float32Col(dflt=UNSET_FLOAT)
    combined_eqw_rest_lya_err = tables.Float32Col(dflt=UNSET_FLOAT)

    spectral_slope = tables.Float32Col(dflt=UNSET_FLOAT)
    spectral_slope_err = tables.Float32Col(dflt=0.0)

    ccd_adjacent_mag = tables.Float32Col(dflt=99.9) #ccd adjacent, single fiber brightest mag
    central_single_fiber_mag = tables.Float32Col(dflt=99.9) #ccd adjacent, single fiber brightest mag
    ffsky_subtraction = tables.BoolCol(dflt=False) #if true, this detection used the full-field sky subtraction

    classification_labels = tables.StringCol(itemsize=32,dflt="")

    #add filter colors
    color_ug = tables.Float32Col(shape=(3,),dflt=[np.nan,np.nan,np.nan] ) #as color, blue max, red_max
    color_ur = tables.Float32Col(shape=(3,),dflt=[np.nan,np.nan,np.nan] ) #as color, blue max, red_max
    color_gr = tables.Float32Col(shape=(3,),dflt=[np.nan,np.nan,np.nan] ) #as color, blue max, red_max

    obs_total_exptime = tables.Float32Col(dflt=UNSET_FLOAT)  #added 0.9.2
    obs_num_dithers = tables.Int8Col(dflt=0)  #added 0.9.2


class SpectraLines(tables.IsDescription):
    detectid = tables.Int64Col(pos=0)  # unique HETDEX detection ID 1e9+
    wavelength = tables.Float32Col(dflt=UNSET_FLOAT,pos=1)
    type = tables.Int32Col(dflt=UNSET_INT,pos=2) # 1 = emission, 0 = unknown, -1 = absorbtion
    score = tables.Float32Col(dflt=UNSET_FLOAT,pos=3)
    sn = tables.Float32Col(dflt=UNSET_FLOAT,pos=4)
    chi2 = tables.Float32Col(dflt=UNSET_FLOAT,pos=5)
    flux_line = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_line_err = tables.Float32Col(dflt=UNSET_FLOAT)
    used = tables.BoolCol(dflt=False) #True if used in the reported multiline solution
    sol_num = tables.Int32Col(dflt=-1) #-1 is unset, 0 is simple scan, no solution, 1+ = solution number in
                                       #decreasing score order

    sigma = tables.Float32Col(dflt=UNSET_FLOAT)
    sigma_err = tables.Float32Col(dflt=UNSET_FLOAT)
    continuum = tables.Float32Col(dflt=UNSET_FLOAT)
    continuum_err = tables.Float32Col(dflt=UNSET_FLOAT)

class Fiber2DCutouts (tables.IsDescription):

    detectid = tables.Int64Col(pos=0)
    ###needs to be array of 4 (top 4 fibers)
    fiber_id = tables.StringCol(shape=(4,),itemsize=38, pos=1)
    distance = StdFloatCol(shape=(4,),pos=2)
    weight = StdFloatCol(shape=(4,),pos=3)
    wavelength = tables.Float32Col(shape=(49,),pos=4) #needs the extra sigfigs
    img_sum = StdFloatCol(shape=(9,49),pos=5)
    img_arr = StdFloatCol(shape=(4, 9, 49),pos=6)
    err_arr = StdFloatCol(shape=(4, 9, 49),pos=7)


#these next two tables decide on 16 or 32 at runtime, so can't use the predefs
class CalibratedSpectra16(tables.IsDescription):
    detectid = tables.Int64Col(pos=0)  # unique HETDEX detection ID 1e9+
    #wavelength = tables.Float32Col(shape=(1036,),pos=1) #skip it
    flux = StdFloatCol(shape=(1036,),pos=2 )  #from Float32Col
    flux_err = StdFloatCol(shape=(1036,),pos=3)  #from Float32Col
    #new 0.9.0
    dust_corr = HalfFloatCol(shape=(1036,), pos=4) #from Float32Col #dust multiplier (normally already applied) to flux and flux_err
    aperture_radius = HalfFloatCol(dflt=UNSET_FLOAT, pos=5) #from Float32Col
    sky_background = tables.Int8Col(dflt=UNSET_INT,pos=6) #from Int32  ... basically 0 for local 1 for ffsky
    num_fibers = tables.Int16Col(dflt=UNSET_INT,pos=7) #Int8 should actually be enough, but just incase a big aperture comes in ....

class CalibratedSpectra32(tables.IsDescription):
    detectid = tables.Int64Col(pos=0)  # unique HETDEX detection ID 1e9+
    #wavelength = tables.Float32Col(shape=(1036,),pos=1) #skip it
    flux = tables.Float32Col(shape=(1036,),pos=2 )  #from Float32Col
    flux_err = tables.Float32Col(shape=(1036,),pos=3)  #from Float32Col
    #new 0.9.0
    dust_corr = tables.Float16Col(shape=(1036,),pos=4) #from Float32Col #dust multiplier (normally already applied) to flux and flux_err
    aperture_radius = tables.Float16Col(dflt=UNSET_FLOAT,pos=5) #from Float32Col
    sky_background = tables.Int8Col(dflt=UNSET_INT,pos=6) #from Int32  ... basically 0 for local 1 for ffsky
    num_fibers = tables.Int16Col(dflt=UNSET_INT,pos=7) #Int8 should actually be enough, but just incase a big aperture comes in ....


class Aperture(tables.IsDescription):
    #one entry per aperture photometry collected
    detectid = tables.Int64Col(pos=0)
    ra = tables.Float32Col(pos=1,dflt=UNSET_FLOAT) #was aperture_ra
    dec = tables.Float32Col(pos=2,dflt=UNSET_FLOAT) #was aperture_dec
    catalog_name = tables.StringCol(itemsize=16)
    filter_name = tables.StringCol(itemsize=16)
    image_depth_mag = HalfFloatCol(dflt=99.9)
    pixel_scale = StdFloatCol(dflt=UNSET_FLOAT)
    radius = StdFloatCol(dflt=UNSET_FLOAT) #in arcsec , #was aperture_radius
    mag = HalfFloatCol(dflt=UNSET_FLOAT) #was aperture_mag
    mag_err = HalfFloatCol(dflt=UNSET_FLOAT) #was  aperture_mag_err
    mag_dered = HalfFloatCol(dflt=UNSET_FLOAT) # #added 0.9.2
    aperture_area_pix = tables.Float32Col(dflt=UNSET_FLOAT) #pixels  #leave at 32
    sky_area_pix = tables.Float32Col(dflt=UNSET_FLOAT) #pixels  #leave at 32
    eqw_rest_lya = StdFloatCol(dflt=UNSET_FLOAT) #was  aperture_eqw_rest_lya
    eqw_rest_lya_err = StdFloatCol(dflt=UNSET_FLOAT) #was  aperture_eqw_rest_lya_err
    plae = StdFloatCol(dflt=UNSET_NAN) #was  aperture_plae
    plae_max = StdFloatCol(dflt=UNSET_NAN) #was  aperture_plae_max
    plae_min = StdFloatCol(dflt=UNSET_NAN) #was  aperture_plae_min
    aperture_cts = tables.Float32Col(dflt=UNSET_FLOAT) #was aperture_counts  #leave at 32
    sky_cts = tables.Float32Col(dflt=UNSET_FLOAT)  #leave at 32
    sky_average = tables.Float32Col(dflt=UNSET_FLOAT)  #leave at 32

class ElixerApertures(tables.IsDescription): #perhaps keep just the selected ones??
    #one entry per aperture photometry collected
    detectid = tables.Int64Col(pos=0)
    ra = tables.Float32Col(pos=1,dflt=UNSET_FLOAT) #decimal degrees of center
    dec = tables.Float32Col(pos=2,dflt=UNSET_FLOAT)
    catalog_name = tables.StringCol(itemsize=16)
    filter_name = tables.StringCol(itemsize=16)
    image_depth_mag = HalfFloatCol(dflt=99.9)
    pixel_scale = StdFloatCol(dflt=UNSET_FLOAT) #arcsec/pixel
    selected = tables.BoolCol(dflt=False) #if True this is the object used for the aperture PLAE/OII, etc (see above table)
    radius = StdFloatCol(dflt=UNSET_FLOAT) #major axis (diameter) 'a' in arcsec
    mag = HalfFloatCol(dflt=UNSET_FLOAT)
    mag_err = HalfFloatCol(dflt=UNSET_FLOAT)
    mag_dered = HalfFloatCol(dflt=UNSET_FLOAT)  #added 0.9.2
    background_cts = tables.Float32Col(dflt=UNSET_FLOAT) #sky_average
    background_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flags = tables.Int32Col(dflt=0) #aperture flags
    image_flags = tables.Int64Col(dflt=0) #separate from the aperture flags, these are ties to the image reduction pipeline


class ExtractedObjects(tables.IsDescription): #keep all, I think
    #one entry per aperture photometry collected
    detectid = tables.Int64Col(pos=0)
    ra = tables.Float32Col(pos=1,dflt=UNSET_FLOAT) #decimal degrees of center
    dec = tables.Float32Col(pos=2,dflt=UNSET_FLOAT)
    catalog_name = tables.StringCol(itemsize=16)
    filter_name = tables.StringCol(itemsize=16)
    image_depth_mag = HalfFloatCol(dflt=99.9)
    pixel_scale = StdFloatCol(dflt=UNSET_FLOAT) #arcsec/pixel
    selected = tables.BoolCol(dflt=False) #if True this is the object used for the aperture PLAE/OII, etc (see above table)
    major = StdFloatCol(dflt=UNSET_FLOAT) #major axis (diameter) 'a' in arcsec
    minor = StdFloatCol(dflt=UNSET_FLOAT) #'b'
    theta = StdFloatCol(dflt=0.0) #radians counter-clockwise from x-axis
    mag = HalfFloatCol(dflt=UNSET_FLOAT)
    mag_err = HalfFloatCol(dflt=UNSET_FLOAT)
    mag_dered = HalfFloatCol(dflt=UNSET_FLOAT)  #added 0.9.2

    background_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    background_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flags = tables.Int32Col(dflt=0)
    dist_curve = StdFloatCol(dflt=UNSET_FLOAT)
    dist_baryctr = StdFloatCol(dflt=UNSET_FLOAT)
    image_flags = tables.Int64Col(dflt=0) #separate from the aperture flags, these are ties to the image reduction pipeline

    fixed_aper_radius = StdFloatCol(dflt=UNSET_FLOAT)
    fixed_aper_mag = StdFloatCol(dflt=UNSET_FLOAT)
    fixed_aper_mag_err = StdFloatCol(dflt=UNSET_FLOAT)
    fixed_aper_mag_dered = StdFloatCol(dflt=UNSET_FLOAT)  #added 0.9.2

    fixed_aper_flux_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    fixed_aper_flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    fixed_aper_flags = tables.Int32Col(dflt=0)

    # flag ... bit mask
    # 01 sep.OBJ_MERGED	      object is result of deblending
    # 02 sep.OBJ_TRUNC	      object is truncated at image boundary
    # 08 sep.OBJ_SINGU	      x, y fully correlated in object
    # 10 sep.APER_TRUNC	      aperture truncated at image boundary
    # 20 sep.APER_HASMASKED	  aperture contains one or more masked pixels
    # 40 sep.APER_ALLMASKED	  aperture contains only masked pixels
    # 80 sep.APER_NONPOSITIVE aperture sum is negative in kron_radius



class CatalogMatch(tables.IsDescription):
    # one entry per catalog bid target
    detectid = tables.Int64Col(pos=0)
    ra = tables.Float32Col(pos=1,dflt=UNSET_FLOAT) #was cat_ra
    dec = tables.Float32Col(pos=2,dflt=UNSET_FLOAT) #was cat_dec
    selected = tables.BoolCol(pos=3,dflt=False)
    catalog_name = tables.StringCol(itemsize=16)
    filter_name = tables.StringCol(itemsize=16)
    match_num = tables.Int16Col(dflt=-1)
    separation = StdFloatCol(dflt=UNSET_FLOAT) #in arcsec
    prob_match = HalfFloatCol(dflt=UNSET_FLOAT) #in arcsec
    specz = HalfFloatCol(dflt=UNSET_FLOAT) #was cat_specz
    photz = HalfFloatCol(dflt=UNSET_FLOAT) #was cat_photz
    flux = tables.Float32Col(dflt=UNSET_FLOAT) #was  cat_flux
    flux_err = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_flux_err
    mag = HalfFloatCol(dflt=UNSET_FLOAT) #was  cat_mag
    mag_err = HalfFloatCol(dflt=UNSET_FLOAT) #was cat_mag_err
    eqw_rest_lya = StdFloatCol(dflt=UNSET_FLOAT) #was cat_eqw_rest_lya
    eqw_rest_lya_err = StdFloatCol(dflt=UNSET_FLOAT) #was cat_eqw_rest_lya_err
    plae = StdFloatCol(dflt=UNSET_FLOAT) #was  cat_plae
    plae_max = StdFloatCol(dflt=UNSET_FLOAT) #was cat_plae_max
    plae_min = StdFloatCol(dflt=UNSET_FLOAT) #was cat_plae_min

    #maybe add in the PDF of the photz ... not sure how big
    #to make the columns ... needs to be fixed, but might
    #vary catalog to catalog
    #cat_photz_pdf_z = tables.Float32Col(1036, )
    #cat_photz_pdf_p = tables.Float32Col(1036, )



class VIRUSShot(tables.IsDescription): #Shot table
    """
    cloned from HETDEX_API create_shot_hdf5.py
    """
    shotid = tables.Int64Col(pos=0)
    date = tables.Int32Col(pos=3)
    obsid = tables.Int32Col(pos=4)
    objid = tables.StringCol((18), pos=2)
    field = tables.StringCol((12), pos=1)
    ra = tables.Float64Col(pos=5)
    dec = tables.Float64Col(pos=6)
    pa = tables.Float64Col(pos=7)
    n_ifu = tables.Int32Col(pos=8)
    datevobs = tables.StringCol((12))
    trajcra = tables.Float32Col()
    trajcdec = tables.Float32Col()
    trajcpa = tables.Float32Col()
    structaz = tables.Float32Col()
    time = tables.StringCol(7)
    ambtemp = tables.Float32Col()
    humidity = tables.Float32Col()
    dewpoint = tables.Float32Col()
    pressure = tables.Float32Col()
    expnum = tables.Int32Col((3))
    exptime = tables.Float32Col((3))
    darktime = tables.Float32Col((3))
    mjd = tables.Float32Col((3))
    fwhm_virus = tables.Float32Col(pos=9)
    fwhm_virus_err = tables.Float32Col(pos=10)
    nstars_fit_fwhm = tables.Int32Col()
    # relflux_guider = tb.Float32Col((3),pos=13)
    relflux_virus = tables.Float32Col((3), pos=14)
    response_4540 = tables.Float32Col(pos=11)  # normalized for 360s
    xditherpos = tables.Float32Col((3))
    yditherpos = tables.Float32Col((3))
    xoffset = tables.Float32Col((3))
    yoffset = tables.Float32Col((3))
    xrms = tables.Float32Col((3))
    yrms = tables.Float32Col((3))
    nstars_fit = tables.Int32Col((3))

    obsind = tables.Int32Col()


#these next two tables decide on 16 or 32 at runtime, so can't use the predefs
class VIRUSFiber16(tables.IsDescription): #uses Float16 where possibly
    """
    cloned from HETDEX_API create_shot_hdf5.py

    remove redundant columns (for this purpose ... single shot only)
    reduce the unneccesary precision to save space

    """
    #obsind = tables.Int32Col() #usuall only 3 digits, but I suppose technically could be up to 9, so leave as is?
                                #BUT this is redundant with the shot table since this is just one shot
                                #also redundant with the multiframe string

    fiber_id = tables.StringCol((38), pos=0)
    #flag = tables.Int32Col(dflt=1, pos=1)  #1 is "good" (copied from VIRUSFiberIndexWithFlags aka /FiberIndex
    ra = tables.Float32Col(pos=1) #may still want this precision, float16 really only gives 4 decimals here
    dec = tables.Float32Col(pos=2)
    multiframe = tables.StringCol((20),pos=3) #this is redundant with fiber_id (which contains the multiframe)

    fibnum = tables.Int8Col(pos=4) # only runs 1 to 112
    fibidx = tables.Int8Col(pos=5)  # this the index on the amp (e.g. 0 to 111) #really redundant with fibnum-1
    healpix = tables.Int64Col(pos=6) # from VIRUSFiberIndex


    ifux = tables.Float32Col()
    ifuy = tables.Float32Col()
    fpx = tables.Float32Col()
    fpy = tables.Float32Col()


    calfib = StdFloatCol((1036,))
    calfibe = StdFloatCol((1036,))
    calfib_ffsky = StdFloatCol((1036,)) #could consider dropping this tone too

    ifuslot = tables.StringCol(3)
    ifuid = tables.StringCol(3)
    specid = tables.StringCol(3)
    contid = tables.StringCol(8)
    amp = tables.StringCol(2)
    expnum = tables.Int32Col()

    #cloned from VIRUSFiberIndex
    flag = tables.Int8Col(dflt=1) #was Int32 #pos=19) #1 is "good" * note: 1 IFF all others are 1, but is not (currently) a bitmap
    flag_badamp = tables.Int8Col(dflt=1) #pos=20) #1 is "good"
    flag_badfib = tables.Int8Col(dflt=1)#pos=21) #1 is "good"
    flag_meteor = tables.Int8Col(dflt=1)#pos=22) #1 is "good"
    flag_satellite = tables.Int8Col(dflt=1)#pos=23) #1 is "good"
    flag_largegal = tables.Int8Col(dflt=1)#pos=24) #1 is "good"
    flag_shot = tables.Int8Col(dflt=1) #pos=25) #1 is "good"
    flag_throughput = tables.Int8Col(dflt=1) #pos=26) #1 is "good"

    #consider removing these to save longterm storage
    #if needed, would re-run
    spectrum = tables.Float16Col((1032,))
    wavelength = tables.Float32Col((1032,))

    fiber_to_fiber = tables.Float16Col((1032,)) #re-extraction in HETDEX_API needs it

    #
    chi2 = tables.Float16Col((1032,))
    rms = tables.Float16Col((1032,))
    # calfib_counts = tables.Float16ColCol((1036,))
    # calfibe_counts = tables.Float16ColCol((1036,))
    #
    trace = tables.Float16Col((1032,))
    sky_subtracted = tables.Float16Col((1032,))
    sky_spectrum = tables.Float16Col((1032,))
    #error1D = tables.Float16ColCol((1032,))

class VIRUSFiber32(tables.IsDescription): #same as VIRUSFiber but uses Float32 instead of Float16 if needed
        """
        cloned from HETDEX_API create_shot_hdf5.py

        remove redundant columns (for this purpose ... single shot only)
        reduce the unneccesary precision to save space

        """
        # obsind = tables.Int32Col() #usuall only 3 digits, but I suppose technically could be up to 9, so leave as is?
        # BUT this is redundant with the shot table since this is just one shot
        # also redundant with the multiframe string

        fiber_id = tables.StringCol((38), pos=0)
        # flag = tables.Int32Col(dflt=1, pos=1)  #1 is "good" (copied from VIRUSFiberIndexWithFlags aka /FiberIndex
        ra = tables.Float32Col(pos=1)  # may still want this precision, float16 really only gives 4 decimals here
        dec = tables.Float32Col(pos=2)
        multiframe = tables.StringCol((20), pos=3)  # this is redundant with fiber_id (which contains the multiframe)

        fibnum = tables.Int8Col(pos=4)  # only runs 1 to 112
        fibidx = tables.Int8Col(pos=5)  # this the index on the amp (e.g. 0 to 111) #really redundant with fibnum-1
        healpix = tables.Int64Col(pos=6)  # from VIRUSFiberIndex

        ifux = tables.Float32Col()
        ifuy = tables.Float32Col()
        fpx = tables.Float32Col()
        fpy = tables.Float32Col()

        calfib = tables.Float32Col((1036,))
        calfibe = tables.Float32Col((1036,))
        calfib_ffsky = tables.Float32Col((1036,))  # could consider dropping this tone too

        ifuslot = tables.StringCol(3)
        ifuid = tables.StringCol(3)
        specid = tables.StringCol(3)
        contid = tables.StringCol(8)
        amp = tables.StringCol(2)
        expnum = tables.Int32Col()

        # cloned from VIRUSFiberIndex
        flag = tables.Int8Col(dflt=1)  # was Int32 #pos=19) #1 is "good" * note: 1 IFF all others are 1, but is not (currently) a bitmap
        flag_badamp = tables.Int8Col(dflt=1)  # pos=20) #1 is "good"
        flag_badfib = tables.Int8Col(dflt=1)  # pos=21) #1 is "good"
        flag_meteor = tables.Int8Col(dflt=1)  # pos=22) #1 is "good"
        flag_satellite = tables.Int8Col(dflt=1)  # pos=23) #1 is "good"
        flag_largegal = tables.Int8Col(dflt=1)  # pos=24) #1 is "good"
        flag_shot = tables.Int8Col(dflt=1)  # pos=25) #1 is "good"
        flag_throughput = tables.Int8Col(dflt=1)  # pos=26) #1 is "good"

        # consider removing these to save longterm storage
        # if needed, would re-run
        spectrum = tables.Float32Col((1032,))
        wavelength = tables.Float32Col((1032,))

        fiber_to_fiber = tables.Float32Col((1032,)) #re-extraction in HETDEX_API needs it as part of the masking

        #
        chi2 = tables.Float32Col((1032,))
        rms = tables.Float32Col((1032,))
        # calfib_counts = tables.Float16Col((1036,))
        # calfibe_counts = tables.Float16Col((1036,))
        #
        trace = tables.Float32Col((1032,))
        sky_subtracted = tables.Float32Col((1032,))
        sky_spectrum = tables.Float32Col((1032,))
        # error1D = tables.Float16ColCol((1032,))


#VIRUSImage, very expensive ... more than half total shot.h5 is wrapped up here. Will skip.
#limit only to critical ones and limit to 16 bit?
class VIRUSImage16(tables.IsDescription):
    obsind = tables.Int32Col()
    multiframe = tables.StringCol((20), pos=0)
    image = tables.Float16Col((1032, 1032))  #elixer uses
    error = tables.Float16Col((1032, 1032)) #elixer uses
    clean_image = tables.Float16Col((1032, 1032)) #elixer uses ??
    ifuslot = tables.StringCol(3)
    ifuid = tables.StringCol(3)
    specid = tables.StringCol(3)
    contid = tables.StringCol(8)
    amp = tables.StringCol(2)
    expnum = tables.Int16Col()

class VIRUSImage32(tables.IsDescription):
    obsind = tables.Int32Col()
    multiframe = tables.StringCol((20), pos=0)
    image = tables.Float32Col((1032, 1032))
    error = tables.Float32Col((1032, 1032))
    clean_image = tables.Float32Col((1032, 1032))
    ifuslot = tables.StringCol(3)
    ifuid = tables.StringCol(3)
    specid = tables.StringCol(3)
    contid = tables.StringCol(8)
    amp = tables.StringCol(2)
    expnum = tables.Int16Col()

class AmpStats(tables.IsDescription):
    """
    cloned from HETDEX_AI ampstats.py

    keep only strictly necessary info
    """
    # shotid ... do not need shotid
    multiframe = tables.StringCol(itemsize=20, pos=0)
    expnum = tables.Int8Col(pos=1)  #
    status = tables.Int32Col(pos=2)  # a status indicator, TBD ... could be a value or a bitmapped mask (-1 bad, 0 unchecked, 1 good?)
    flag = tables.Int32Col(pos=3)

    #the below are the base information used to set the status
    #could re-run the reduction to recreate and keep these out, but there are no arrays here and < 1000 amps
    # so is not that much data ... just keep it?
    im_median = tables.Float32Col()
    mask_fraction = tables.Float32Col()
    avg = tables.Float32Col()
    scale = tables.Float32Col()
    chi2fib_med = tables.Float32Col()
    frac_c2 = tables.Float32Col()
    frac_0 = tables.Float32Col()
    n_lo = tables.Int32Col()
    avg_orig = tables.Float32Col()
    sky_sub_rms = tables.Float32Col()
    sky_sub_rms_rel = tables.Float32Col()
    sky_sub_rms_median = tables.Float32Col()
    dither_relflux = tables.Float32Col()
    norm = tables.Float32Col()
    kchi = tables.Float32Col()
    n_cont = tables.Int32Col()



class CalfibDQ(tables.IsDescription):
    """
    cloned from HETDEX_API create_fiber_mask_hdf5.py
    """

    fiber_id = tables.StringCol((38), pos=0)
    calfib_dq = tables.Int16Col((1036,),pos=1)

# for reference, these are the bitmapped mask values
# class CALFIB_DQ(BitFlagNameMap):
#     MAIN = 1                  # 0x0001
#     FTF = 2                   # 0x0002
#     CHI2FIB = 4               # 0x0004
#     BADPIX = 8                # 0x0008
#     BADAMP = 16               # 0x0010
#     LARGEGAL = 32             # 0x0020
#     METEOR = 64               # 0x0040
#     BADSHOT = 128             # 0x0080
#     THROUGHPUT = 256          # 0x0100
#     BADFIB = 512              # 0x0200
#     SATELLITE = 1024          # 0x0400
#     BADCAL = 2048             # 0x0800
#     PIXMASK = 4096            # 0x1000

# now combined with Fibers
class VIRUSFiberIndexWithFlags(tables.IsDescription):
    """
    cloned from HETDEX_API create_fiber_index_hdr5.py

    remove redundant info to save space

    note: this is the alternate form of VIRUSFiberIndex
          it is used in the same way but (in original HETDEX) adds the flag information
          It both forms are used to create a group/table named FiberIndex, with the Single Shot Reduction
          using the VIRUSFiberIndexWithFlags variant instead

    *** NOTE: there will likely be some code modifications needed in HETDEX_API to look for the flag in the new place
              if it is not found in the old group/table

    """
    #minimum necessary
    multiframe = tables.StringCol((20), pos=0) #redundant with Fibers table
    ra = tables.Float32Col(pos=1) #redundant with Fibers table
    dec = tables.Float32Col(pos=2) #redundant with Fibers table
    fiber_id = tables.StringCol((38),pos=3)
    healpix = tables.Int64Col(pos=4)

    #do we need these? maybe for speedy selection
    amp = tables.StringCol(2,pos=5) #redundant with Fibers table
    date = tables.Int64Col(pos=6) #redundant with Fibers table
    datevobs = tables.StringCol((12),pos=7) #redundant with Fibers table
    expnum = tables.Int32Col(pos=8) #redundant with Fibers table
    fibidx = tables.Int8Col() #redundant with Fibers table
    fibnum = tables.Int32Col(pos=10) #redundant with Fibers table
    fpx = tables.Float32Col(pos=11) #redundant with Fibers table
    fpy = tables.Float32Col(pos=12) #redundant with Fibers table
    ifuslot = tables.StringCol(3,pos=13) #redundant with Fibers table
    ifuid = tables.StringCol(3,pos=14) #redundant with Fibers table
    ifux = tables.Float32Col(pos=15) #redundant with Fibers table
    ifuy = tables.Float32Col(pos=16) #redundant with Fibers table
    shotid = tables.Int64Col(pos=17) #redundant with Fibers table
    specid = tables.StringCol(3,pos=18) #redundant with Fibers table

    flag = tables.Int32Col(dflt=1) #was Int32 #pos=19) #1 is "good" * note: 1 IFF all others are 1, but is not (currently) a bitmap
    flag_badamp = tables.Int8Col(dflt=1) #pos=20) #1 is "good"
    flag_badfib = tables.Int8Col(dflt=1)#pos=21) #1 is "good"
    flag_meteor = tables.Int8Col(dflt=1)#pos=22) #1 is "good"
    flag_satellite = tables.Int8Col(dflt=1)#pos=23) #1 is "good"
    flag_largegal = tables.Int8Col(dflt=1)#pos=24) #1 is "good"
    flag_shot = tables.Int8Col(dflt=1) #pos=25) #1 is "good"
    flag_throughput = tables.Int8Col(dflt=1) #pos=26) #1 is "good"


class ThroughputTable (tables.IsDescription):

    wavelength = tables.Int32Col(pos=0)
    throughput = tables.Int32Col(pos=1)
    tp_low = tables.Int32Col(pos=2)
    tp_high = tables.Int32Col(pos=3)
    rat_poly = tables.Int32Col(pos=4)
    tp_gband = tables.Int32Col(pos=5)


class QualityAssessment(tables.IsDescription):
    expnum = tables.Int32Col(pos=0)
    xoffset = tables.Float32Col(pos=1)
    yoffset = tables.Float32Col(pos=2)
    xrms = tables.Float32Col(pos=3)
    yrms = tables.Float32Col(pos=4)
    nstars = tables.Float32Col(pos=5)


class NominalVals(tables.IsDescription):
    expnum = tables.Int32Col(pos=0)
    ra = tables.Float32Col(pos=1)
    dec = tables.Float32Col(pos=2)
    parangle = tables.Float32Col(pos=3)
    x_dither_pos = tables.Float32Col(pos=4)
    y_dither_pos = tables.Float32Col(pos=5)
    relflux_virus = tables.Float32Col(pos=6)
    als_filename = tables.StringCol((23))


class Fplane(tables.IsDescription):
    ifuslot = tables.Int32Col(pos=0)
    fpx = tables.Float32Col(pos=1)
    fpy = tables.Float32Col(pos=2)
    specid = tables.Int32Col(pos=3)
    specslot = tables.Int32Col(pos=4)
    ifuid = tables.Int32Col(pos=5)
    ifurot = tables.Float32Col(pos=6)
    platesc = tables.Float32Col(pos=7)

class Dithall(tables.IsDescription):
    ra = tables.Float32Col(pos=0)
    dec = tables.Float32Col(pos=1)
    ifuslot = tables.StringCol((3))
    XS = tables.Float32Col()
    YS = tables.Float32Col()
    xfplane = tables.Float32Col()
    yfplane = tables.Float32Col()
    multifits = tables.StringCol((28))
    timestamp = tables.StringCol((17))
    exposure = tables.StringCol((5))

class PositionOffsets(tables.IsDescription):
    xoffset = tables.Float32Col(pos=0)
    yoffset = tables.Float32Col(pos=1)
    ra_dex = tables.Float32Col(pos=2)
    dec_dex = tables.Float32Col(pos=3)
    ra_cat = tables.Float32Col(pos=4)
    dec_cat = tables.Float32Col(pos=5)
    ifuslot = tables.Int32Col(pos=6)


class StarCatalog(tables.IsDescription):
    ignore = tables.Int64Col(pos=0)
    star_ID = tables.StringCol((28), pos=1)
    ra_cat = tables.Float64Col(pos=2)
    dec_cat = tables.Float64Col(pos=3)
    u = tables.Float64Col(pos=4)
    g = tables.Float64Col(pos=5)
    r = tables.Float64Col(pos=6)
    i = tables.Float64Col(pos=7)
    z = tables.Float64Col(pos=8)


class StarCatalogMatches(tables.IsDescription):
    RA_det = tables.Float32Col(pos=0)
    DEC_det = tables.Float32Col((28), pos=1)
    IFUSLOT_det = tables.Int32Col(pos=2)
    xifu_det = tables.Float32Col(pos=3)
    yifu_det = tables.Float32Col(pos=4)
    ifuslot_det = tables.Int32Col(pos=5)
    RA_cat = tables.Float32Col(pos=6)
    DEC_cat = tables.Float32Col(pos=7)
    IFUSLOT_cat = tables.Int32Col(pos=8)
    xifu_cat = tables.Float32Col(pos=9)
    yifu_cat = tables.Float32Col(pos=10)
    ifuslot_cat = tables.Int32Col(pos=11)

def build_ssr_shot_h5(shot_fn, elixer_fn=None):#, outfn=None):
    """
    combine and trim the original <shot.h5> and <elixer.h5> files
    if elixer_fn is None, just trim the <shot.h5>


    :param shot_fn: MUST provide
    :param elixer_fn: can be None

    ## no, want the naming fixed, so do not offer to the user, (can always change in the file system later)
    ##:param outfn: this is the output file name, should provide but defaults to ssr+<shot_fn>
    :return:
    """

    global FATAL_EXIT, elixer_h5_force

    fileh = None
    shot_h5 = None
    elixer_h5 = None
    datevshot = None
    outfn = None

    try:
        datevshot = os.path.basename(shot_fn).replace(".h5", "")
    except:
        pass

    try:
        #check that the one or two input files exist and can be opened
        if shot_fn is not None:
            try:
                shot_h5 = tables.open_file(shot_fn,mode='r')
            except:
                log.critical(f"[{datevshot}] Failed to open input shot.h5 file: {shot_fn}",exc_info=True)
                return

        if shot_h5 is None: #this MUST always be provided
            log.critical(f"[{datevshot}] Fatal. Cannot open shot.h5 file: {shot_fn}")
            return

        if elixer_fn is not None:
            try:
                elixer_h5 = tables.open_file(elixer_fn,mode='r')
            except:
                log.critical(f"[{datevshot}] Fatal. Failed to open input elixer.h5 file: {elixer_fn}",exc_info=True)
                return

            if elixer_h5 is None:
                log.critical(f"[{datevshot}] Fatal. Failed to open input elixer.h5 file: {elixer_fn}", exc_info=True)
                return

            #check supported version
            elixer_h5_version = elixer_h5.root.Version.read(field='version')[0]
            if elixer_h5_version not in SUPPORTED_ELIXER_H5_VERSIONS:
                if elixer_h5_force:
                    log.critical(f"[{datevshot}] WARNING! Input elixer h5 version {elixer_h5_version} not supported."
                                 f" Forcing attempt to build ssr h5 anyway. {elixer_fn} ")
                else:
                    log.critical(f"[{datevshot}] Fatal. Input elixer h5 version {elixer_h5_version} not supported. {elixer_fn} ")
                    FATAL_EXIT = 1
                    return


        outfn = "ssr_" + os.path.basename(shot_fn)

        log.debug(f"[{datevshot}] ***START*** Creating new SingleShot Reduction HDF5 catalog (%s)" % (outfn))

        fileh = tables.open_file(outfn, 'w', 'SingleShot Reduction Catalog')

        vtb = fileh.create_table(fileh.root, 'Version', Version,
                                 'Version Table')

        row = vtb.row
        row['version'] = __version__
        row['version_pytables'] = tables.__version__
        row.append()
        vtb.flush()

        #######################################
        # Shot Table
        #######################################

        fileh.create_table(fileh.root, 'Shot', VIRUSShot, 'Shot Summary Table')

        #just one row for Shot table, so just do a direct copy
        fileh.root.Shot.append(shot_h5.root.Shot.read())
        fileh.root.Shot.flush()

        ############################
        # Calibration/Throughput
        # very small table with some limited use
        # (including for compatibility with analysis notebook)
        ############################

        try:
            group = fileh.create_group(fileh.root, "Calibration", "HETDEX Calibration Info")
            groupThroughput = fileh.create_group(group, "Throughput", "Throughput Curve")
            fileh.create_table(fileh.root.Calibration.Throughput, "throughput", ThroughputTable, 'Throughput Curve')

            for row in shot_h5.root.Calibration.Throughput.throughput.read():
                # for row in shot_h5.root.Data.Fibers.read():
                new_row = fileh.root.Calibration.Throughput.throughput.row
                # go over some columns individual since want to change some types
                # directy copy columns:
                for col in fileh.root.Calibration.Throughput.throughput.colnames:
                    new_row[col] = row[col].astype(np.float32) #they are float64, needlessly

                new_row.append()

            fileh.root.Calibration.Throughput.throughput.flush()


            #and the image /Calibration/Throughput/tp_png (ImageArray(1684, 1190, 4))
            try:
                if "tp_png" in [node.name for node in shot_h5.root.Calibration.Throughput]:
                    pngimarr = shot_h5.root.Calibration.Throughput.tp_png.read()
                    pngim = fileh.create_array(groupThroughput, "tp_png", pngimarr,filter=COMPRESSION_FILTER)
                    pngim.attrs["CLASS"] = "IMAGE"
                    pngim.flush()
                else:
                    print(f"[{datevshot}] Non fatal. Will ignore. tp_png not found in root.Calibration.Throughput")

            except:
                print(f"[{datevshot}] Non fatal. Will ignore. Fail on Calibration/Throughput/tp_png: {traceback.format_exc()}")

        except:
            #log.debug("Non fatal. Will ignore. Fail on Calibration/Throughput/throughput table", exc_info=True)
            print(f"[{datevshot}] Non fatal. Will ignore. Fail on Calibration/Throughput: {traceback.format_exc()}")


        #######################################
        # FullSkyModel (1 per exposure)
        #######################################
        #notice: these are different from the others are are separate groups, each
        # with a single 2D array
        # We do, however, want to reduce from float64 to float32
        # the array is wave (which needs float32) and counts (which could be float16 and be okay)
        #  but we can't mix types this way unless we change the format

        print(f"[{datevshot}] Importing FullSkyModels ... ", flush=True)

        groupFullSkyModel = fileh.create_group(fileh.root, 'FullSkyModel', 'FullSkyModel')
        #node_names = [node.name for node in shot_h5.root.FullSkyModel._f_list_nodes()]
        for node in shot_h5.root.FullSkyModel:
            n = shot_h5.get_node(node)
            sky_array = n.read().astype(np.float32)
            fileh.create_array(groupFullSkyModel,n.name,sky_array)

        fileh.flush()



        #######################################
        # Fibers
        #######################################

        print(f"[{datevshot}] Importing Fiber data ... ", flush=True)

        #put under "Data" for some compatibility with pathing in HETDEX_API
        _ = fileh.create_group(fileh.root, "Data", "VIRUS Fiber Data") # was assigned to group_data previously

        if StdFloatCol == FullFloatCol:
            use32 = True
        else:
            use32 = False
            #in order of most likely to need float32
            f16_check_fields = ["calfib","calfib_ffsky","spectrum","sky_subtracted","sky_spectrum","trace"]
            print(f"[{datevshot}] Checking for float16 compatibility ... ")
            for fd in f16_check_fields:
                mx = np.max(shot_h5.root.Data.Fibers.read(field=fd))
                if mx > threshold_16:
                    print(f"[{datevshot}] {fd} (at least) requires float32: mx {mx}")
                    use32 = True
                    #print(f"{fd} : min {np.min(shot_h5.root.Data.Fibers.read(field=fd))} , max {np.max(shot_h5.root.Data.Fibers.read(field=fd))}")
                    break #for a test, don't break ... want to see them all

            f16_fields_to_clip = ["chi2"]

            # use32 = False
            # #mx = np.abs(shot_h5.root.Data.Fibers.read(field="calfib")).max()
            # mx = np.max(shot_h5.root.Data.Fibers.read(field="calfib"))
            # #mn = np.min(shot_h5.root.Data.Fibers.read(field="calfib"))
            # if mx > threshold_16:
            #     use32 = True
            # else:
            #     #mx = np.abs(shot_h5.root.Data.Fibers.read(field="calfib_ffsky")).max()
            #     mx = np.max(shot_h5.root.Data.Fibers.read(field="calfib_ffsky"))
            #     if mx > threshold_16:
            #         use32 = True
            #     else:
            #         mx = np.max(shot_h5.root.Data.Fibers.read(field="calfibe"))  # can only be positive
            #         if mx > threshold_16:
            #             use32 = True
            #         else:
            #             mx = np.max(shot_h5.root.Data.Fibers.read(field="fiber_to_fiber"))
            #             if mx > threshold_16:
            #                 use32 = True

        #note: turning on fitlers=COMPRESSION_FITLER can save about 500MB per file BUT
        #      reading the fibers is 1.5 to 3x longer and fibers are what are hit most
        if use32:
            print(f"[{datevshot}] Using Float32 for VIRUSFibers",flush=True)
            fileh.create_table(fileh.root, 'Fibers', VIRUSFiber32, 'Fiber Summary Table')#,filters=COMPRESSION_FILTER)
        else:
            print(f"[{datevshot}] Using Float16 for VIRUSFibers",flush=True)
            fileh.create_table(fileh.root, 'Fibers', VIRUSFiber16, 'Fiber Summary Table')#,filters=COMPRESSION_FILTER)

        #this will also be softlinked to root.Data.Fibers
        #fileh.create_table(group_data, 'Fibers', VIRUSFiber, 'Fiber Summary Table')
       # fileh.create_table(fileh.root, 'Fibers', VIRUSFiber16, 'Fiber Summary Table')


        #need the flags
        # Combine with FiberIndex and then make a softlink
        fi_dict = dict(zip(shot_h5.root.FiberIndex.read(field="fiber_id"), shot_h5.root.FiberIndex.read()))

        # to instead keep the FiberIndex table, even with its redundant info
        # flags = shot_h5.root.FiberIndex.read(field="flag")
        # fiber_ids = list(shot_h5.root.FiberIndex.read(field="fiber_id"))
        # flag_dict = dict(zip(fiber_ids,flags))
        # del flags
        # del fiber_ids

        #many for fibers, so iterate
        for row in tqdm(shot_h5.root.Data.Fibers.read(),disable=not SHOW_TQDM):
        #for row in shot_h5.root.Data.Fibers.read():
            new_row = fileh.root.Fibers.row
            #go over some columns individual since want to change some types
            #directy copy columns:
            for col in ['multiframe','fiber_id','ifux','ifuy','fpx','fpy','ra','dec',
                        'ifuslot','ifuid','specid','contid','amp','fiber_to_fiber']:
                new_row[col] = row[col]

            #special treatment, mostly float32 to float16
            #expnum (int32 to int16)
            #fibnum (int32 to int8)
            #calfib (float32 array to float16 array) #Float16Col((1036,))
            #calfibe (float32 array to float16 array) #Float16Col((1036,))
            #calfib_ffsky (float32 array to float16 array) #Float16Col((1036,))

            new_row['fibnum'] = row['fibnum'].astype(np.int8) #1 to 112
            new_row['fibidx'] = row['fibidx'].astype(np.int8) #just 0 to 111, the index on the amp
            new_row['expnum'] = row['expnum'].astype(np.int16) #Maybe we'd have more than 15 exposures (prob not, though)
            if use32:
                new_row['calfib'] = row['calfib'].astype(np.float32)
                new_row['calfibe'] = row['calfibe'].astype(np.float32)
                new_row['calfib_ffsky'] = row['calfib_ffsky'].astype(np.float32)

                #questionable ones
                new_row['wavelength'] = row['wavelength'].astype(np.float32)  # always 32-bit
                new_row['spectrum'] = row['spectrum'].astype(np.float32)
                new_row['chi2'] = row['chi2'].astype(np.float32)
                new_row['rms'] = row['rms'].astype(np.float32)
                new_row['trace'] = row['trace'].astype(np.float32)
                new_row['sky_subtracted'] = row['sky_subtracted'].astype(np.float32)
                new_row['sky_spectrum'] = row['sky_spectrum'].astype(np.float32)
            else:
                # could have stupidly negative values (always errorneous and since already checked for max okay, apply a clip to be safe)
                #new_row['calfib'] = row['calfib'].astype(StdFloat)
                new_row['calfib'] = np.clip(row['calfib'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                new_row['calfibe'] = np.clip(row['calfibe'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                #could have stupidly negative values (always errorneous and since already checked for max okay, apply a clip to be safe)
                #new_row['calfib_ffsky'] = row['calfib_ffsky'].astype(StdFloat)
                new_row['calfib_ffsky'] = np.clip(row['calfib_ffsky'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)

                # questionable ones
                new_row['wavelength'] = row['wavelength'].astype(np.float32) #always 32-bit
                new_row['spectrum'] = np.clip(row['spectrum'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                new_row['chi2'] = np.clip(row['chi2'],a_min=-1*clip_float16,a_max=clip_float16).astype(StdFloat)
                new_row['rms'] = np.clip(row['rms'],a_min=-1*clip_float16,a_max=clip_float16).astype(StdFloat)
                new_row['trace'] = np.clip(row['trace'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                new_row['sky_subtracted'] = np.clip(row['sky_subtracted'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                new_row['sky_spectrum'] = np.clip(row['sky_spectrum'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)


            #add in the corresponding FiberIndex
            fi_row = fi_dict[new_row['fiber_id']]

            for col in ['healpix','flag','flag_badamp','flag_badfib','flag_meteor','flag_satellite',
                        'flag_largegal','flag_shot','flag_throughput']:
                new_row[col] = fi_row[col]


            # #need the flags ... NO, see above note about FiberIndex table
            # try:
            #     new_row['flag'] = flag_dict[new_row['fiber_id']]
            # except:
            #     log.warning("Could not retrieve fibers flag", exc_info=True)
            #     new_row['flag'] = 1

            new_row.append()

        fileh.root.Fibers.flush()

        #set the index
        try:
            fileh.root.Fibers.cols.fiber_id.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on fibers table:", exc_info=True)
        try:
            fileh.root.Fibers.cols.multiframe.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on fibers table:", exc_info=True)
        try:
            fileh.root.Fibers.cols.ra.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on fibers table:", exc_info=True)
        try:
            fileh.root.Fibers.cols.dec.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on fibers table:", exc_info=True)
        try:
            fileh.root.Fibers.cols.healpix.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on fibers table:",exc_info=True)

        fileh.root.Fibers.flush()

        #print("Trying softlink ....")
        shot_h5.create_soft_link(fileh.root.Data, 'Fibers', target=fileh.root.Fibers)



        #####################################
        # FiberIndex  --
        # DO NOT combine with Fibers Table ... the way HETDEX_API reads in (all at once) is too costly
        #     with the other columns of the Fibers Table
        ######################################

        fileh.create_table(fileh.root, 'FiberIndex', VIRUSFiberIndexWithFlags, 'FiberIndex Table')
        #note, below we will create a softlink under root.Data.FiberIndex for compatibility

        # many for fibers, so iterate
        # root.Data.FiberIndex or root.FiberIndex  (root.Data.FiberIndex Does not have the flags)
        print(f"[{datevshot}] Importing Fiber data ... ", flush=True)
        for row in tqdm(shot_h5.root.FiberIndex.read(),disable=not SHOW_TQDM):
            # for row in shot_h5.root.Data.Fibers.read():
            new_row = fileh.root.FiberIndex.row
            # go over some columns individual since want to change some types
            # directy copy columns:
            for col in fileh.root.FiberIndex.colnames:
                #['fiber_id','healpix','flag_badamp','flag_badfib','flag_meteor','flag_satellite',
                #        'flag_largegal','flag_shot','flag_throughput']:
                new_row[col] = row[col]

            #changed the size of this one
            #new_row['flag'] = row['flag'].astype(np.int8)

            new_row.append()

        fileh.root.FiberIndex.flush()

        # set the index
        try:
            fileh.root.FiberIndex.cols.fiber_id.create_csindex()
            fileh.root.FiberIndex.cols.healpix.create_csindex()
            #fileh.root.FiberIndex.cols.ra.create_csindex() #HETDEX_API is not using ra, dec as index
            #fileh.root.FiberIndex.cols.dec.create_csindex()
        except:
            log.debug("Index fail on FiberIndex table", exc_info=True)

        fileh.root.FiberIndex.flush()

        #create a softlink for compatibility ????
        print(f"[{datevshot}] Trying softlink to FiberIndex....")
        shot_h5.create_soft_link(fileh.root.Data, 'FiberIndex', target=fileh.root.FiberIndex)





        ##################################
        # VIRUSImages
        # ... VERY expensive (more than half the total size of the shot.h5
        #      will skip. Does not impact re-extraction, BUT you cannot build elixer reports without them
        #                    and without other data.
        ##################################
        if not minimum_h5:
            print(f"[{datevshot}] Importing VIRUS CCD image data ... ", flush=True)

            if StdFloatCol == FullFloatCol:
                use32 = True
                clip_error = 0.
            else:
                use32 = False
                clip_error = 0.
                #note: should not expect grossly negative values, but may be a scaling issue?
                #for now, at least, allow absolute value check
                mx = np.abs(shot_h5.root.Data.Images.read(field="image")).max()
                if mx > image_threshold_16:
                    use32 = True
                    log.debug(f"[{datevshot}] Data.Images.image requires float32: mx {mx}")
                else:
                    mx = np.abs(shot_h5.root.Data.Images.read(field="clean_image")).max()
                    if mx > image_threshold_16:
                        use32 = True
                        log.debug(f"[{datevshot}] Data.Images.clean_image requires float32: mx {mx}")
                    else:
                        mx = np.max(shot_h5.root.Data.Images.read(field="error"))  # can only be positive
                        if mx > image_threshold_16:
                            #want to set max value to 65504.00
                            #use32 = True
                            #log.debug(f"Data.Images.error requires float32: mx {mx}")
                            clip_error = 65500.00 #max is 655504
                            log.debug(f"[{datevshot}] Data.Images.error clipping to {clip_error} : mx {mx}")

            if use32:
                print(f"[{datevshot}] Using Float32 for VIRUSImages")
                fileh.create_table(fileh.root.Data, 'Images', VIRUSImage32, 'VIRUS CCD Image Data',filters=COMPRESSION_FILTER)
            else:
                print(f"[{datevshot}] Using Float16 for VIRUSImages")
                fileh.create_table(fileh.root.Data, 'Images', VIRUSImage16, 'VIRUS CCD Image Data',filters=COMPRESSION_FILTER)

            for row in tqdm(shot_h5.root.Data.Images.read(),disable=not SHOW_TQDM):
                # for row in shot_h5.root.Data.Fibers.read():
                new_row = fileh.root.Data.Images.row
                # go over some columns individual since want to change some types
                # directy copy columns:
                for col in ['obsind', 'multiframe', 'ifuslot', 'ifuid', 'specid', 'contid', 'amp']:
                    new_row[col] = row[col]

                # special treatment, mostly float32 to float16
                new_row['expnum'] = row['expnum'].astype(np.int16)  # Maybe we'd have more than 15 exposures (prob not, though)
                if use32:
                    new_row['image'] = row['image'].astype(np.float32)
                    new_row['error'] = row['error'].astype(np.float32)
                    new_row['clean_image'] = row['clean_image'].astype(np.float32)
                else:
                    new_row['image'] = row['image'].astype(np.float16)
                    #new_row['error'] = row['error'].astype(np.float16)
                    if clip_error > 0:
                        new_row['error'] = np.clip(row['error'],a_min=-1*clip_error,a_max=clip_error).astype(np.float16)
                    new_row['clean_image'] = row['clean_image'].astype(np.float16)


                new_row.append()

            fileh.root.Data.Images.flush()

            # set the index
            try:
                fileh.root.Data.Images.cols.multiframe.create_csindex()
            except:
                log.debug(f"Index fail on Data.Images table", exc_info=True)

            fileh.root.Data.Images.flush()

        #####################################
        # CalfibDQ
        ######################################

        fileh.create_table(fileh.root, 'CalfibDQ', CalfibDQ, 'Fiber per-Wavelength Flags Table')
        print(f"[{datevshot}] Importing CalfibDQ (Fiber per-wavelength flags) data ...",flush=True)
        copy_cols = fileh.root.CalfibDQ.colnames

        for row in tqdm(shot_h5.root.CalfibDQ.read(),disable=not SHOW_TQDM):
            new_row = fileh.root.CalfibDQ.row
            # go over some columns individual since want to change some types
            # directy copy columns:
            for col in copy_cols:
                new_row[col] = row[col]

            new_row.append()

        fileh.root.CalfibDQ.flush()
        # set the index
        try:
            fileh.root.CalfibDQ.cols.fiber_id.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on CalfibDQ table", exc_info=True)

        fileh.root.CalfibDQ.flush()



        #####################################
        # AmpStats
        ######################################

        fileh.create_table(fileh.root, 'AmpStats', AmpStats, 'Amp Stats Table')
        print(f"[{datevshot}] Importing AmpStats data ...",flush=True)
        copy_cols = fileh.root.AmpStats.colnames
        copy_cols.remove('expnum')

        for row in tqdm(shot_h5.root.AmpStats.read(),disable=not SHOW_TQDM):
            new_row = fileh.root.AmpStats.row
            # go over some columns individual since want to change some types
            # directy copy columns:
            for col in copy_cols:
                new_row[col] = row[col]

            # changed the size of this one
            new_row['expnum'] = row['expnum'].astype(np.int8)

            new_row.append()

        fileh.root.AmpStats.flush()
        # set the index
        try:
            fileh.root.AmpStats.cols.multiframe.create_csindex()
        except:
            log.debug(f"[{datevshot}] Index fail on AmpStats table", exc_info=True)

        fileh.root.AmpStats.flush()



        ######################################
        # Astrometry
        # Images
        # NOTE: unlike the elixer images, these are 1 per leaf table
        ######################################

        if full_h5 or Include_CoaddImage_matched or Include_CoaddImage_ifugrid or Include_CoaddImage_skyview:
            groupAstrometry = fileh.create_group(fileh.root, 'Astrometry', 'Astrometry Info')

            print(f"[{datevshot}] Importing Collapsed Images ...", flush=True)
            try:
                #create the group with compression
                groupCoadd = fileh.create_group(groupAstrometry, 'CoaddImages', 'Coadd Images',filters=COMPRESSION_FILTER)

                for node in shot_h5.root.Astrometry.CoaddImages:
                    n = shot_h5.get_node(node)

                    if (n.name[0:3] == "mat" and Include_CoaddImage_matched) or \
                       (n.name[0:3] == "png" and Include_CoaddImage_ifugrid) or \
                       (n.name[0:3] == "exp" and Include_CoaddImage_skyview):
                        try:
                            img_array = n.read().astype(np.float32)
                            #compression is automtic due to group creation with filters defined
                            fileh.create_array(groupCoadd, n.name, img_array)
                        except:
                            log.debug(f"[{datevshot}] Exception adding {n}", exc_info=True)

                fileh.flush()
            except:
                log.debug(f"[{datevshot}] Exception adding Collapsed Image(s)", exc_info=True)



            if full_h5: #the rest of the Astrometry

                try:
                    print(f"[{datevshot}] Astrometry.NominalVals")
                    fileh.create_table(groupAstrometry, 'NominalVals', NominalVals, 'Nominal Values')
                    copy_cols = fileh.root.Astrometry.NominalVals.colnames
                    for row in tqdm(shot_h5.root.Astrometry.NominalVals.read()):#, disable=not SHOW_TQDM):
                        new_row = fileh.root.Astrometry.NominalVals.row
                        # go over some columns individual since want to change some types
                        # directy copy columns:
                        for col in copy_cols:
                            new_row[col] = row[col]

                        new_row.append()

                    fileh.root.Astrometry.NominalVals.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.NominalVals", exc_info=True)

                try:
                    print(f"[{datevshot}] Astrometry.QA")
                    fileh.create_table(groupAstrometry, 'QA', QualityAssessment, 'Quality Assessment')
                    copy_cols = fileh.root.Astrometry.QA.colnames
                    for row in tqdm(shot_h5.root.Astrometry.QA.read()):#, disable=not SHOW_TQDM):
                        new_row = fileh.root.Astrometry.QA.row
                        # go over some columns individual since want to change some types
                        # directy copy columns:
                        for col in copy_cols:
                            new_row[col] = row[col]

                        new_row.append()

                    fileh.root.Astrometry.QA.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.QA", exc_info=True)

                #print(f"[{datevshot}] Deliberately skipping Astrometry.ShuffleCfg. This is by design.")

                try:
                    blob = shot_h5.root.Astrometry.ShuffleCfg.read()
                    if len(blob) > 0:
                        fileh.create_array(groupAstrometry, 'ShuffleCfg', blob)
                        shot_h5.root.Astrometry.ShuffleCfg.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.ShuffleCfg", exc_info=True)


                try:
                    print(f"[{datevshot}] Astrometry.StarCatalog")
                    fileh.create_table(groupAstrometry, 'StarCatalog', StarCatalog, 'StarCatalog')
                    copy_cols = fileh.root.Astrometry.StarCatalog.colnames
                    for row in tqdm(shot_h5.root.Astrometry.StarCatalog.read()):#, disable=not SHOW_TQDM):
                        new_row = fileh.root.Astrometry.StarCatalog.row
                        # go over some columns individual since want to change some types
                        # directy copy columns:
                        for col in copy_cols:
                            new_row[col] = row[col]

                        new_row.append()

                    fileh.root.Astrometry.StarCatalog.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.StarCatalog", exc_info=True)


                try:
                    print(f"[{datevshot}] Astrometry.fplane")
                    fileh.create_table(groupAstrometry, 'fplane', Fplane, 'fplane')
                    int_cols = ['ifuslot','specid','specslot','ifuid']
                    float_cols = ['fpx','fpy','ifurot','platesc']
                    for row in tqdm(shot_h5.root.Astrometry.fplane.read()):#, disable=not SHOW_TQDM):
                        new_row = fileh.root.Astrometry.fplane.row
                        # go over some columns individual since want to change some types
                        # directy copy columns:
                        for col in int_cols:
                            new_row[col] = row[col].astype(np.int32)
                        for col in float_cols:
                            new_row[col] = row[col].astype(np.float32)

                        new_row.append()

                    fileh.root.Astrometry.fplane.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.fplane", exc_info=True)


                #CatalogMatches ... depends on the number of exposures
                try:
                    print(f"[{datevshot}] Astrometry.CatalogMatches.exp??")
                    subgroup = fileh.create_group(fileh.root.Astrometry, 'CatalogMatches', 'Match Catalog Info')

                    #notice: there is redundant info IFUSLOT_det vs ifuslot_det
                    #                            and IFUSLOT_cat vs ifuslot_cat
                    # ASSUME this is for some compatibility? overall size impact is negligible so just leave as is
                    int_cols = ['IFUSLOT_det','ifuslot_det','IFUSLOT_cat','ifuslot_cat']
                    float_cols = ['RA_det','DEC_det','xifu_det','yifu_det','RA_cat','DEC_cat','xifu_cat','yifu_cat']

                    for node in shot_h5.root.Astrometry.CatalogMatches: #this is "exp01", "exp02", ...
                        scmtb = fileh.create_table(subgroup, node.name, StarCatalogMatches, 'Match Catalog Info')
                        for row in node:
                            new_row = scmtb.row

                            #there is no need for the int64 and float64, but int16 and float16 is not enough
                            #so cast to 32bit
                            for col in int_cols:
                                new_row[col] = np.int32(row[col])#.astype(np.int32)
                            for col in float_cols:
                                new_row[col] = np.float32(row[col])#.astype(np.float32)

                        scmtb.flush()

                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.CatalogMatches.exp??", exc_info=True)

                # Dithall ... depends on the number of exposures
                try:
                    print(f"[{datevshot}] Astrometry.Dithall.exp??")
                    subgroup = fileh.create_group(fileh.root.Astrometry, 'Dithall', 'Fiber Astrometry Info')

                    for node in shot_h5.root.Astrometry.Dithall:  # this is "exp01", "exp02", ...
                        datb = fileh.create_table(subgroup, node.name, Dithall, 'Fiber Astrometry Info')
                        for row in node:
                            new_row = datb.row
                            copy_cols = node.colnames

                            for col in copy_cols:
                                new_row[col] = row[col]

                        datb.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.Dithall.exp??", exc_info=True)


                # PositionOffsets ... depends on the number of exposures
                try:
                    print(f"[{datevshot}] Astrometry.PositionOffsets.exp??")
                    subgroup = fileh.create_group(fileh.root.Astrometry, 'PositionOffsets', 'Offset in star matches')
                    int_cols = ['ifuslot']
                    float_cols = ['xoffset','yoffset','ra_dex','dec_dex','ra_cat','dec_cat']

                    for node in shot_h5.root.Astrometry.PositionOffsets:  # this is "exp01", "exp02", ...
                        potb = fileh.create_table(subgroup, node.name, PositionOffsets, 'Offset in star matches')
                        for row in node:
                            new_row = potb.row

                            for col in float_cols:
                                new_row[col] = np.float32(row[col])#.astype(np.float32)
                            for col in int_cols:
                                new_row[col] = np.int32(row[col])#.astype(np.int32)

                        potb.flush()
                except:
                    log.debug(f"[{datevshot}] Exception adding Astrometry.PositionOffsets.exp??", exc_info=True)

        #done with the shot.h5 file now
        shot_h5.close()








        if elixer_h5 is not None:
            print(f"[{datevshot}] Importing ELiXer data ... ",flush=True)

            fileh.create_table(fileh.root, 'Detections', Detections,
                               'Detection Summary Table')

            fileh.create_table(fileh.root, 'SpectraLines', SpectraLines,
                               'Identified SpectraLines Table')

            #since there is a choice about 16 or 32 bit, this is moved later
            # fileh.create_table(fileh.root, 'CalibratedSpectra', CalibratedSpectra,
            #                    'PSF Weighted Spectra Table')

            fileh.create_table(fileh.root, 'Aperture', Aperture,
                               'Aperture Photometry Table') # mostly a g and r aperture, sometimes more

            fileh.create_table(fileh.root, 'CatalogMatch', CatalogMatch,
                               'Archival Catalog Matched Objected Table')

            fileh.create_table(fileh.root, 'ExtractedObjects', ExtractedObjects,
                               'SourceExtractor Objects Table')  # multiple filters, many objects

            fileh.create_table(fileh.root, 'ElixerApertures', ElixerApertures,
                               'Circular Apertures Table')  # mostly a g and r aperture, sometimes more

            fileh.create_table(fileh.root, 'Fiber2DCutouts', Fiber2DCutouts,
                               '2D CCD Cutouts Around Emission Line Table')

            ###########################################################################
            #Detections all direct copy but not all columns
            ###########################################################################
            skip_cols = np.array(['shotid', 'obsid', 'seeing_fwhm', 'response', 'fieldname'])
            keep_cols = np.setdiff1d(elixer_h5.root.Detections.colnames, skip_cols)
            for row in elixer_h5.root.Detections.read():
                new_row = fileh.root.Detections.row

                # go over some columns individual since want to change some types
                # directy copy columns:
                for col in keep_cols:
                    new_row[col] = row[col]
                new_row.append()

            fileh.root.Detections.flush()
            #set the index
            try:
                fileh.root.Detections.cols.detectid.create_csindex()
                fileh.root.Detections.cols.ra.create_csindex()
                fileh.root.Detections.cols.dec.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on Detections table", exc_info=True)

            fileh.root.Detections.flush()



            ###########################################################################
            #Fiber2DCutouts all direct copy with some resolution change
            ###########################################################################
            try:
                keep_cols = elixer_h5.root.Fiber2DCutouts.colnames

                for row in elixer_h5.root.Fiber2DCutouts.read():
                    new_row = fileh.root.Fiber2DCutouts.row

                    for col in keep_cols:
                        new_row[col] = row[col]
                    new_row.append()

                fileh.root.Fiber2DCutouts.flush()
                #set the index
                try:
                    fileh.root.Fiber2DCutouts.cols.detectid.create_csindex()
                except:
                    log.debug(f"[{datevshot}] Index fail on Fiber2DCutouts table", exc_info=True)

                fileh.root.Fiber2DCutouts.flush()
            except:
                log.info(f"[{datevshot}] Fiber2DCutouts table not found.",  exc_info=True)

            #############################################################################
            #SpectraLines ... unchanged (note: could reduce some precision, but there
            #                 are relatively few entries so not worth much storage
            #############################################################################
            fileh.root.SpectraLines.append(elixer_h5.root.SpectraLines.read())
            fileh.root.SpectraLines.flush()
            #set the index
            try:
                fileh.root.SpectraLines.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on SpectraLines table", exc_info=True)

            fileh.root.SpectraLines.flush()

            #############################################################################
            # CalibratedSpectra
            #############################################################################

            #need to check for maximum size of flux or flux_err

            if StdFloatCol == FullFloatCol:
                use32 = True
            else:
                use32 = False
                #mx = np.abs(elixer_h5.root.CalibratedSpectra.read(field="flux")).max()
                mx = np.max(elixer_h5.root.CalibratedSpectra.read(field="flux"))
                if mx > threshold_16:
                    use32 = True
                else:
                    mx = np.max(elixer_h5.root.CalibratedSpectra.read(field="flux_err")) #can only be positive
                    if mx > threshold_16:
                        use32 = True

            if use32:
                print(f"[{datevshot}] Using float32 for CalibratedSpectra",flush=True)
                fileh.create_table(fileh.root, 'CalibratedSpectra', CalibratedSpectra32,
                               'PSF Weighted Spectra Table')
            else:
                print(f"[{datevshot}] Using float16 for CalibratedSpectra",flush=True)
                fileh.create_table(fileh.root, 'CalibratedSpectra', CalibratedSpectra16,
                               'PSF Weighted Spectra Table')

            for row in elixer_h5.root.CalibratedSpectra.read():
                new_row = fileh.root.CalibratedSpectra.row

                new_row['detectid'] = row['detectid']
                # skip wavelength entirely and set the float32 to float16
                if use32:
                    new_row['flux'] = row['flux'].astype(np.float32)
                    new_row['flux_err'] = row['flux_err'].astype(np.float32)
                else:
                    # could have stupidly negative values (always errorneous and since already checked for max okay, apply a clip to be safe)
                    #new_row['flux'] = row['flux'].astype(StdFloat)
                    new_row['flux'] = np.clip(row['flux'], a_min=-1 * clip_float16, a_max=clip_float16).astype(StdFloat)
                    new_row['flux_err'] = row['flux_err'].astype(StdFloat)

                new_row['dust_corr'] = row['dust_corr'].astype(StdFloat)

                #aperture needs help too (not -10000)
                if row['aperture_radius'] < 0:
                    new_row['aperture_radius'] = -1.0
                else:
                    new_row['aperture_radius'] = row['aperture_radius'].astype(StdFloat)
                if row['sky_background'] < 0:
                    new_row['sky_background'] = -1
                else:
                    new_row['sky_background'] = 1 if row['sky_background'] > 0 else 0 #really a boolean
                new_row['num_fibers'] = row['num_fibers'].astype(np.int16)

                new_row.append()

            fileh.root.CalibratedSpectra.flush()
            #set the index
            try:
                fileh.root.CalibratedSpectra.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on CalibratedSpectra table", exc_info=True)

            fileh.root.CalibratedSpectra.flush()


            ###############################################
            # Aperture
            #################################################
            for row in elixer_h5.root.Aperture.read():
                new_row = fileh.root.Aperture.row

                new_row['detectid'] = row['detectid']
                new_row['ra'] = row['ra']
                new_row['dec'] = row['dec']
                new_row['catalog_name'] = row['catalog_name']
                new_row['filter_name'] = row['filter_name']
                new_row['image_depth_mag'] = row['image_depth_mag'].astype(HalfFloat)
                new_row['pixel_scale'] = row['pixel_scale'].astype(StdFloat)
                new_row['radius'] = row['radius'].astype(StdFloat)
                new_row['mag'] = row['mag'].astype(HalfFloat)
                new_row['mag_err'] = row['mag_err'].astype(HalfFloat)
                try:
                    new_row['mag_dered'] = row['mag_dered'].astype(HalfFloat)
                except:
                    new_row['mag_dered'] = 99.9

                new_row['aperture_area_pix'] = row['aperture_area_pix'] #.astype(np.float16) #could be very many pix
                new_row['sky_area_pix'] = row['sky_area_pix'] #.astype(np.float16) #could be very many pix
                new_row['eqw_rest_lya'] = row['eqw_rest_lya'].astype(StdFloat)
                new_row['eqw_rest_lya_err'] = row['eqw_rest_lya_err'].astype(StdFloat)
                new_row['plae'] = row['plae'].astype(StdFloat)
                new_row['plae_max'] = row['plae_max'].astype(StdFloat)
                new_row['plae_min'] = row['plae_min'].astype(StdFloat)

                new_row['aperture_cts'] = row['aperture_cts'] #could be very many counts
                new_row['sky_cts'] = row['sky_cts']   #could be very many counts
                new_row['sky_average'] = row['sky_average'] #could be very many counts

                new_row.append()

            fileh.root.Aperture.flush()
            #set the index
            try:
                fileh.root.Aperture.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on Aperture table", exc_info=True)

            fileh.root.Aperture.flush()



            ##########################################
            # ElixerApertures
            # Only keep the one selected aperture size ?
            ##########################################
            for row in elixer_h5.root.ElixerApertures.read_where("selected==True"):
                new_row = fileh.root.ElixerApertures.row

                new_row['detectid'] = row['detectid']
                new_row['ra'] = row['ra']
                new_row['dec'] = row['dec']
                new_row['catalog_name'] = row['catalog_name']
                new_row['filter_name'] = row['filter_name']
                new_row['image_depth_mag'] = row['image_depth_mag'].astype(HalfFloat)
                new_row['pixel_scale'] = row['pixel_scale'].astype(StdFloat)
                new_row['selected'] = row['selected']
                new_row['radius'] = row['radius'].astype(StdFloat)
                new_row['mag'] = row['mag'].astype(HalfFloat)
                new_row['mag_err'] = row['mag_err'].astype(HalfFloat)
                try:
                    new_row['mag_dered'] = row['mag_dered'].astype(HalfFloat)
                except:
                    new_row['mag_dered'] = 99.9

                new_row['background_cts'] = row['background_cts'] #.astype(np.float16) #could be very many pix
                new_row['background_err'] = row['background_err'] #.astype(np.float16) #could be very many pix


                new_row['flux_cts'] = row['flux_cts'] #could be very many counts
                new_row['flux_err'] = row['flux_err']   #could be very many counts
                new_row['flags'] = row['flags']
                new_row['image_flags'] = row['image_flags']
                new_row.append()

            fileh.root.ElixerApertures.flush()
            #set the index
            try:
                fileh.root.ElixerApertures.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on ElixerApertures table", exc_info=True)

            fileh.root.ElixerApertures.flush()



            ##################################################
            # ExtractedObjects
            ##################################################
            for row in elixer_h5.root.ExtractedObjects.read_where("selected==True"):
                new_row = fileh.root.ExtractedObjects.row

                new_row['detectid'] = row['detectid']
                new_row['ra'] = row['ra']
                new_row['dec'] = row['dec']
                new_row['catalog_name'] = row['catalog_name']
                new_row['filter_name'] = row['filter_name']
                new_row['image_depth_mag'] = row['image_depth_mag'].astype(HalfFloat)
                new_row['pixel_scale'] = row['pixel_scale'].astype(StdFloat)
                new_row['selected'] = row['selected']

                new_row['major'] = row['major'].astype(StdFloat)
                new_row['minor'] = row['minor'].astype(StdFloat)
                new_row['theta'] = row['theta'].astype(StdFloat)
                new_row['mag'] = row['mag'].astype(HalfFloat)
                new_row['mag_err'] = row['mag_err'].astype(HalfFloat)
                try:
                    new_row['mag_dered'] = row['mag_dered'].astype(HalfFloat)
                except:
                    new_row['mag_dered'] = -99.9

                new_row['background_cts'] = row['background_cts'] #.astype(np.float16) #could be very many pix
                new_row['background_err'] = row['background_err'] #.astype(np.float16) #could be very many pix
                new_row['flux_cts'] = row['flux_cts'] #could be very many counts
                new_row['flux_err'] = row['flux_err']   #could be very many counts
                new_row['flags'] = row['flags']

                new_row['dist_curve'] = row['dist_curve'].astype(StdFloat)
                new_row['dist_baryctr'] = row['dist_baryctr'].astype(StdFloat)

                new_row['image_flags'] = row['image_flags']

                new_row['fixed_aper_radius'] = row['fixed_aper_radius'].astype(StdFloat)
                new_row['fixed_aper_mag'] = row['fixed_aper_mag'].astype(HalfFloat)
                new_row['fixed_aper_mag_err'] = row['fixed_aper_mag_err'].astype(HalfFloat)
                try:
                    new_row['fixed_aper_mag_dered'] = row['fixed_aper_mag_dered'].astype(HalfFloat)
                except:
                    new_row['fixed_aper_mag_dered'] = 99.9

                new_row['fixed_aper_flux_cts'] = row['fixed_aper_flux_cts']
                new_row['fixed_aper_flux_err'] = row['fixed_aper_flux_err']
                new_row['fixed_aper_flags'] = row['fixed_aper_flags']

                new_row.append()

            fileh.root.ExtractedObjects.flush()
            #set the index
            try:
                fileh.root.ExtractedObjects.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on ExtractedObjects table", exc_info=True)

            fileh.root.ExtractedObjects.flush()



            #########################################
            # CatalogMatch
            ########################################
            for row in elixer_h5.root.CatalogMatch.read_where("selected==True"):
                new_row = fileh.root.CatalogMatch.row

                new_row['detectid'] = row['detectid']
                new_row['ra'] = row['ra']
                new_row['dec'] = row['dec']
                new_row['selected'] = row['selected']
                new_row['catalog_name'] = row['catalog_name']
                new_row['filter_name'] = row['filter_name']

                new_row['match_num'] = row['match_num'].astype(np.int16)
                new_row['separation'] = row['separation'].astype(StdFloat)
                new_row['prob_match'] = row['prob_match'].astype(HalfFloat)
                new_row['specz'] = row['specz'].astype(HalfFloat)
                new_row['photz'] = row['photz'].astype(HalfFloat)

                new_row['flux'] = row['flux']
                new_row['flux_err'] = row['flux_err']

                new_row['mag'] = row['mag'].astype(HalfFloat)
                new_row['mag_err'] = row['mag_err'].astype(HalfFloat)
                new_row['eqw_rest_lya'] = row['eqw_rest_lya'].astype(StdFloat)
                new_row['eqw_rest_lya_err'] = row['eqw_rest_lya_err'].astype(StdFloat)
                new_row['plae'] = row['plae'].astype(StdFloat)
                new_row['plae_max'] = row['plae_max'].astype(StdFloat)
                new_row['plae_min'] = row['plae_min'].astype(StdFloat)

                new_row.append()

            fileh.root.CatalogMatch.flush()
            #set the index
            try:
                fileh.root.CatalogMatch.cols.detectid.create_csindex()
            except:
                log.debug(f"[{datevshot}] Index fail on CatalogMatch table", exc_info=True)

            fileh.root.CatalogMatch.flush()



            # done with the elixer.h5 file now
            elixer_h5.close()

    except:
        print(f"[{datevshot}] Exception building new h5 file. {traceback.format_exc()}")



    if fileh is not None:
        fileh.close()

    return outfn
#end build_ssr_shot_h5


def get_max_image(image_path,datevshot="???"):
    """

    :param image_path:
    :return: 3-tuple of (max) shape, np.array of unique 1st dimensions
    """
    # only first dimension is allowed to change
    max1 = 0
    max2 = 0
    max3 = 0

    t1 = []
    # t2 = []
    # t3 = []

    try:
        print(f"[{datevshot}] Checking image sizes ...",flush=True)
        image_fns = sorted(glob.glob(image_path))
        for img_path in tqdm(image_fns,disable=not SHOW_TQDM):
            x1,x2,x3  = np.array(Image.open(img_path)).shape

            max1 = max(max1, x1)
            max2 = max(max2, x2)
            max3 = max(max3, x3)

            t1.append(x1)
            # t2.append(x2)
            # t3.append(x3)

    except:
        print(f"[{datevshot}] Exception: {traceback.format_exc()}")

    # print(f"x1: {np.unique(t1)}")
    # print(f"x2: {np.unique(t2)}")
    # print(f"x3: {np.unique(t3)}")

    unique_d1, unique_ct = np.unique(t1,return_counts=True)

    return (max1,max2,max3), unique_d1, unique_ct




def import_images_earray (shot_h5fn,image_path,group_name,earray_name="image_data"):
    """

    :param shot_h5fn: new ssr shot h5 path and filename
    :param image_path: path to images, include wildcards as will be used with glob
    :return:
    """
    datevshot = shot_h5fn
    try:

        datevshot = os.path.basename(shot_h5fn).replace(".h5", "").replace("ssr_","")
        print(f"[{datevshot}] Importing images: {image_path} to root.{group_name}.{earray_name}*",flush=True)

        #max_shape, unique_d1, unique_ct = get_max_image(image_path,datevshot)

        max_shape, img_dict = get_image_dict(image_path,datevshot)
        unique_d1 = img_dict.keys()
        if unique_d1 is None or len(unique_d1) < 1:
            print(f"[{datevshot}] No matching report images found.")
            return
        unique_ct = np.array([len(img_dict[k]) for k in img_dict.keys()])

        #sort
        unique_d1, unique_ct = zip(*sorted(zip(unique_d1, unique_ct))[::-1]) #decreasing size (this is slightly better on storage)
        #unique_d1, unique_ct = zip(*sorted(zip(unique_d1, unique_ct))) #increasing size

        image_fns = []
        for k in unique_d1: #in same D1 size order
            image_fns += img_dict[k]


        with tables.open_file(shot_h5fn,mode="r+") as h5:  #so will auto close regardless of exit
            dtb = h5.root.Detections

            #image_fns = sorted(glob.glob(image_path)) #instead use order from img_dict above

            # Create a group to store the images #might already exist
            try:
                img_group = h5.create_group(h5.root, group_name)
            except:
                img_group = h5.root.elixer_reports

            img = Image.open(image_fns[0])
            img_dtype = np.array(img).dtype

            # Define atom and filters for compression (e.g., zlib, blosc)
            atom = tables.Atom.from_dtype(img_dtype)
            # Using 'blosc' for efficient lossless compression is common 0 = none, 1=minimum up to 9=maximum compression
            #filters = tables.Filters(complevel=9, complib='blosc')
            #filters = tables.Filters(complevel=9, complib='bzip2', bitshuffle=False, shuffle=False)
            # filters = None #now a global COMPRESSION_FILTER

            # Create a resizable eArray to store the images
            # The shape is (0, ...) to start empty, and the maxshape is (None, ...)
            # to allow the first dimension to grow indefinitely.

            all_ea = []
            all_ea_idx = []
            for d1,img_ct in zip(unique_d1,unique_ct):
                name = earray_name + "_" + str(d1)
                img_shape = (d1,max_shape[1],max_shape[2])
                try:
                    image_array = h5.create_earray(img_group, name, atom,
                                                     shape=(0,) + img_shape,
                                                     expectedrows=img_ct,
                                                     filters=COMPRESSION_FILTER)

                    all_ea.append(image_array)
                    all_ea_idx.append(0)
                except:
                    try:
                        e1 = traceback.format_exc()
                        image_array = h5.get_node(img_group)._f_get_child(name)
                    except:
                        print(f"[{datevshot}] Cannot import images.\nException #1: {e1}\nException #2: {traceback.format_exc()}")

            # Iterate through all image files, resize if needed, and append to the eArray
            total_images = len(image_fns)
            print(f"[{datevshot}] Importing {total_images} images ... ",flush=True)
            for img_path in tqdm(image_fns,disable=not SHOW_TQDM):
                img = Image.open(img_path)
                # Optional: Resize images to a consistent size if necessary
                # img = img.resize((new_width, new_height))
                img_as_array = np.array(img)
                d1 = img_as_array.shape[0]
                i = list(unique_d1).index(d1)
                # if img_as_array.shape[0] < max_shape[0]:
                #     print(f"Resizing {os.path.basename(img_path)} from {img_as_array.shape} to {max_shape}")
                #     img_as_array.resize(max_shape)


                #debug:
                #print(i,d1,img_as_array.shape,all_ea[i].name,all_ea[i]._v_chunkshape,img_path)

                all_ea[i].append([img_as_array])

                #now update Detections
                did = Path(img_path).stem
                nei = False
                if did[-4:]=="_nei":
                    nei = True
                    did = did[:-4] #strip off the _nei
                did = np.int64(did)

                idx = dtb.get_where_list("detectid==did")[0]
                row = dtb.read_where("detectid==did")  # [0]
                if nei:
                    row[0]['h5_neighbor_id'] = d1
                    row[0]['h5_neighbor_idx'] = all_ea_idx[i]
                else:
                    row[0]['h5_report_id'] = d1
                    row[0]['h5_report_idx'] = all_ea_idx[i]
                all_ea_idx[i] += 1
                dtb.modify_rows(start=idx, stop=idx + 1, step=1, rows=row)

            dtb.flush()

            print(f"[{datevshot}] Stored {total_images} images in {shot_h5fn}.")

    except:
        print(f"[{datevshot}] Exception in add_report_images(): {traceback.format_exc()}")


def get_image_dict(image_path,datevshot="???"):
    """

    :param image_path:
    :return: 3-tuple of (max) shape, and dictionary of image paths keyed by the image 1st Dimension length
    """

    max1 = 0
    max2 = 0
    max3 = 0

    img_dict = {}

    try:
        print(f"[{datevshot}] Checking image sizes ...",flush=True)
        image_fns = sorted(glob.glob(image_path))
        for img_path in tqdm(image_fns,disable=not SHOW_TQDM):
            x1,x2,x3  = np.array(Image.open(img_path)).shape

            max1 = max(max1, x1)
            max2 = max(max2, x2)
            max3 = max(max3, x3)

            if x1 in img_dict.keys():
                img_dict[x1].append(img_path)
            else:
                img_dict[x1] = [img_path]

    except:
        print(f"[{datevshot}] Exception: {traceback.format_exc()}",flush=True)

    return (max1,max2,max3), img_dict


def import_images_carray(shot_h5fn,image_path,group_name,carray_name="image_data"):
    """

    :param shot_h5fn: new ssr shot h5 path and filename
    :param image_path: path to images, include wildcards as will be used with glob
    :return:
    """

    datevshot = shot_h5fn
    try:
        datevshot = os.path.basename(shot_h5fn).replace(".h5", "").replace("ssr_", "")
        print(f"[{datevshot}] Importing images: {image_path} to root.{group_name}.{carray_name}*",flush=True)

        max_shape, img_dict = get_image_dict(image_path,datevshot)

        with tables.open_file(shot_h5fn,mode="r+") as h5:  #so will auto close regardless of exit

            dtb = h5.root.Detections

            # Create a group to store the images #might already exist
            try:
                img_group = h5.create_group(h5.root, group_name)
            except:
                img_group = h5.root.elixer_reports


            img = Image.open(img_dict[next(iter(img_dict))][0]) #just need one,so open 1st image in dict
            img_dtype = np.array(img).dtype

            # Define atom and filters for compression (e.g., zlib, blosc)
            atom = tables.Atom.from_dtype(img_dtype)
            # Using 'blosc' for efficient lossless compression is common 0 = none, 1=minimum up to 9=maximum compression
            # better compression net result with shuffle=False (default is True) and bitshuffle=False (default)
            #filters = tables.Filters(complevel=1, complib='bzip2',bitshuffle=False,shuffle=False)
            #filters = None #now a global COMPRESSION_FILTER

            #iterate (in decending order of 1D size) over the img_dict and pre-allocate carrays and then populate
            total_images = sum(len(img_dict[k]) for k in img_dict.keys())
            print(f"[{datevshot}] Importing {total_images} images ... ",flush=True)
            for key in sorted(img_dict.keys())[::-1]:
                name = carray_name + "_" + str(key)
                img_shape = (key, max_shape[1], max_shape[2])
                try:
                    image_array = h5.create_carray(img_group, name, atom,
                                                   shape=(len(img_dict[key]),) + img_shape,
                                                   filters=COMPRESSION_FILTER)
                except:
                    try:
                        e1 = traceback.format_exc()
                        image_array = h5.get_node(img_group)._f_get_child(name)
                    except:
                        print(f"[{datevshot}] Cannot import images.\nException #1: {e1}\nException #2: {traceback.format_exc()}",flush=True)

                #print(f"Importing ID:{key} ...")
                for i, img_path in enumerate(tqdm(img_dict[key],disable=not SHOW_TQDM)):
                    img_as_array = np.array(Image.open(img_path))
                    image_array[i] = [img_as_array]

                    # now update Detections
                    did = Path(img_path).stem
                    nei = False
                    if did[-4:] == "_nei":
                        nei = True
                        did = did[:-4]  # strip off the _nei
                    did = np.int64(did)

                    idx = dtb.get_where_list("detectid==did")[0]
                    row = dtb.read_where("detectid==did")  # [0]
                    if nei:
                        row[0]['h5_neighbor_id'] = key
                        row[0]['h5_neighbor_idx'] = i
                    else:
                        row[0]['h5_report_id'] = key
                        row[0]['h5_report_idx'] = i

                    dtb.modify_rows(start=idx, stop=idx + 1, step=1, rows=row)

            dtb.flush()

            print(f"[{datevshot}] Stored {total_images} images in {shot_h5fn}.")

    except:
        print(f"[{datevshot}] Exception in add_report_images(): {traceback.format_exc()}")

########################################################################
########################################################################
########################################################################
# Main (section)
#   notice: no actual main function
########################################################################
########################################################################
########################################################################

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

print(f"version {__version__}")


if "-help" in args:
    help = """
        Single Shot Reduction HDF5 Aggregation and Compression
        
        usage: python ssr_hdf5.py [switches]
        
        output: new file named ssr_<shot_h5>        
               
        --shot_h5           REQUIRED
            The input shot.h5 file to be compressed (includes the filename). Can be a relative path.

        --bg [#]            optional
            Specify the maximum number of simulatenously running hdf5 constructions
            Note: for a single vm-small node on lonestar6, this should be 3
                  for a single development or normal node, this can be 20-24
                  
        --float32           optional
            Override the compression and force float32 where float16 would otherwise be used due to limited range
            Note: this applies only to float32 data in the original (input) shot.h5 file 
        
        --exclude_ccd       optional
            Do not include the CCD images (3x1032x1032 plus supporting fields for each amp)
            note: in compressed float16 representation this is around 6GB for a typical shot. 
            
        --compression       optional
            Set the compression type and level. Mostly affectes ELiXer report images. The default is (2)
            (1) blosc2 compression at level 9; very fast but lesser compression
                    about 20% faster but 25%+ lesser compression than (2)
            (2) zlib compression at level 1; fast with good compression, 
                   [default] typically around 6 minutes run time per shot 
            (3) bzip2 compression at level 9; very slow but maximum compression
                    about 300% longer to run but 40-50% better compression than (2)
                  
        --elixer_h5         optional
            Path to the <elixer>.h5 file for this shot (includes the filename). Can be a relative path.
            if present, the ELiXer hdf5 tables and the ELiXer report images will be included in the output hdf5 file
            
        --elixer_h5_force   optional
             attempt to force the inclusion of ELiXer data even if the reported elixer version is not explicitly supported
             
        --images            optional
            Path the the ELiXer report (png) images (path only, no filenames). Can be a relative path.
            
    """
    print(help)
    exit(0)

shot_h5_path = None
if "-shot_h5" in args: #path to the shot h5 file
    i = args.index("-shot_h5")
    try:
        shot_h5_path = args[i+1]
    except:
        print(f"Invalid -shot_h5 specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-shot_h5")
else:
    print("Fatal. Must supply --shot_h5 <path to the shot hdf5 file>")

datevshot = os.path.basename(shot_h5_path).replace(".h5","")

elixer_h5_path = None
if "-elixer_h5" in args: #path to the shot h5 file
    i = args.index("-elixer_h5")
    try:
        elixer_h5_path = args[i+1]
    except:
        print(f"Invalid -elixer_h5 specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-elixer_h5")

elixer_h5_force = False
if "-elixer_h5_force" in args: #path to the shot h5 file
    elixer_h5_force = True
    args.remove("-elixer_h5_force")

images_path = None
if "-images" in args:
    i = args.index("-images")
    try:
        images_path = args[i+1]
    except:
        print(f"Invalid -images specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-images")


if "-float32" in args: #force32 bit for fields that were originally 32bit (do not allow down-casting to float16)
    StdFloatCol = FullFloatCol
    args.remove("-float32")

#maybe this is enough ?
# if "-exclude_ccd" in args:
#     exclude_ccd_images = True
#     args.remove("-exclude_ccd")

if "-minimum" in args:
    minimum_h5 = True
    args.remove("-minimum")

if "-bg" in args:
    SHOW_TQDM = False
    i = args.index("-bg")
    try:
        Max_Simultaneous_Shots = int(args[i+1])
    except:
        print(f"Invalid -bg argument specified: {args[i + 1]}")
        exit(-1)
    del args[i + 1]
    args.remove("-bg")
else:
    if not minimum_h5: #cut in half
        Max_Simultaneous_Shots = round(Max_Simultaneous_Shots/2)


#default, level 2
COMPRESSION_FILTER = tables.Filters(complevel=1, complib='zlib', bitshuffle=False, shuffle=False)
if "-compression" in args:
    i = args.index("-compression")
    try:
        compression = int(args[i+1])

        #create the filter to use
        if compression == 1: #least compression (that is still acceptable, but fast)
            #around 5 minutes (300s) for typical ELiXer data ~ 1000 detects and Neighbors
            #around 3.3 GB
            #yes, this one should be at level 9 (basically same time as level1 and a good bit better compression)
            COMPRESSION_FILTER = tables.Filters(complevel=9, complib='blosc2', bitshuffle=False, shuffle=False)
            print(f"[{datevshot}] Using type 1 compression: blosc2 at complvl 9 . Limited compression (lossless), low CPU + time")
        elif compression == 2: #standard
            # around 6.0 minutes (360s) for typical ELiXer data ~ 1000 detects and Neighbors
            # around 2.6 GB
            # leave at complevel 1 (increasing is huge time cost for very little compression improvement)
            COMPRESSION_FILTER = tables.Filters(complevel=1, complib='zlib', bitshuffle=False, shuffle=False)
            print(f"[{datevshot}] Using type 2 compression: zlib at complvl 1 . Good compression, moderate CPU + time")
        elif compression == 3: #maximum
            # around 17 minutes (1000s) for typical ELiXer data ~ 1000 detects and Neighbors
            # around 1.9 GB (1.8Gb at level 9 and 18 minutes)
            # 1 vs 9 is arund 1015s vs 1080s (e.g. +1 minute) for about 5% better compression
            COMPRESSION_FILTER = tables.Filters(complevel=9, complib='bzip2', bitshuffle=False, shuffle=False)
            print(f"[{datevshot}] Using type 3 compression: bzip2 at complvl 9 . Maximum compression, maximum CPU + time")
        else:
            print(f"[{datevshot}] Unexpected --compression value {compression}: Must be in [1,2,3], least compression to most and shortest to longest time cost")
            exit(-1)
    except:
        print(f"[{datevshot}] Invalid -compression specified: {args[i+1]}")
        exit(-1)

    del args[i+1]  # args.pop(0) #remove THIS file
    args.remove("-compression")
else:
    print(f"[{datevshot}] Using [default] type 2 compression: zlib at complvl 1 . Good compression, moderate CPU + time")

if len(args) > 0:
    print(f"[{datevshot}] Unknown remainting args: {args}")



wait_to_run(Max_Simultaneous_Shots, datevshot=datevshot,clean_up=False)

start_time = time.perf_counter()

new_h5_fn = build_ssr_shot_h5(shot_h5_path, elixer_fn=elixer_h5_path)
if images_path is not None and not FATAL_EXIT and new_h5_fn is not None:
    #print("Skipping images import")
    #using the earray and growing
    import_images_earray(new_h5_fn,os.path.join(images_path,"*[0-9].png"),"elixer_reports")
    import_images_earray(new_h5_fn,os.path.join(images_path,"*_nei.png"),"elixer_neighbors")

    #pre-allocating with carray
    #import_images_carray(new_h5_fn,os.path.join(images_path,"*[0-9].png"),"elixer_reports")
    #import_images_carray(new_h5_fn,os.path.join(images_path,"*_nei.png"),"elixer_neighbors")


print(f"[{os.path.basename(shot_h5_path).replace('.h5','')}] Elapsed time: {time.perf_counter() - start_time :0.1f}s")
wait_to_run(Max_Simultaneous_Shots, datevshot=datevshot,clean_up=True)

# example code to fetch images
# main elixer report
# did = np.int64(25112201200666)
# row = dtb.read_where("detectid==did")[0]
# gp = h5.get_node(f"/elixer_reports")
# path = gp._f_get_child(f"image_data_{row['h5_report_id']}")
# idx = row['h5_report_idx']
# Image.fromarray(path[idx]).show()

# neighborhod report
# gp = h5.get_node(f"/elixer_neighbors")
# path = gp._f_get_child(f"image_data_{row['h5_neighbor_id']}")
# idx = row['h5_neighbor_idx']
# Image.fromarray(path[idx]).show()