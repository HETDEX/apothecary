"""
Single Shot Reduction HDF5 file creation

For long term storage of key data for a single shot reduction

Based on the elixer_hdf5.py code AND the HETDEX_API shot h5 code

This file is for a single shot (observation) ONLY. Do NOT commbine shots.

"""


__version__ = '0.1.0'


import numpy as np
import tables
import os
import sys
import time as pytime
import traceback

#for convenience
from elixer import global_config as G

UNSET_FLOAT = -999.999
UNSET_INT = -99999
UNSET_STR = ""
UNSET_NAN = np.nan
SUPPORTED_ELIXER_H5_VERSIONS = [b"0.9.2",]

#logging will just be prints
#this is all single threaded, no real management needed
class PrintLog():
    def __init__(self,logfile=None):
        self.logfile =logfile

    def log(self,msg):
        print(msg,flush=True)
        if self.logfile:
            try:
                with open(self.logfile,"a") as f:
                    f.write(msg+"\n")
            except:
                pass

    def critical(self,msg):
        self.log(msg)

    def error(self,msg):
        self.log(msg)

    def info(self,msg):
        self.log(msg)

    def warning(self,msg):
        self.log(msg)

    def debug(self,msg):
        self.log(msg)

log = PrintLog('h5_merge.log')

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

    elixer_version = tables.StringCol(itemsize=16,pos=2) #version of elixer that generated this detection report
    elixer_datetime = tables.StringCol(itemsize=21,pos=3) #YYYY-MM-DD hh:mm:ss


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
    ra = tables.Float32Col(dflt=UNSET_FLOAT,pos=4)
    dec = tables.Float32Col(dflt=UNSET_FLOAT,pos=5)
    wavelength_obs = tables.Float32Col(dflt=UNSET_FLOAT,pos=6)
    wavelength_obs_err = tables.Float32Col(dflt=UNSET_FLOAT,pos=7)
    apcor_4500 = tables.Float32Col(dflt=UNSET_FLOAT)

    z_best = tables.Float32Col(dflt=-1.0,pos=8)
    z_best_pz = tables.Float32Col(dflt=0.0,pos=9)

    z_best_plya_thresh = tables.Float32Col(dflt=-1.0, pos=10)
    z_best_2 = tables.Float32Col(dflt=-1.0, pos=11)
    z_best_pz_2 = tables.Float32Col(dflt=0.0, pos=12)
    z_best_plya_thresh_2 = tables.Float32Col(dflt=-1.0, pos=13)
    z_best_3 = tables.Float32Col(dflt=-1.0, pos=14)
    z_best_pz_3 = tables.Float32Col(dflt=0.0, pos=15)
    z_best_plya_thresh_3 = tables.Float32Col(dflt=-1.0, pos=16)

    flags = tables.Int32Col(dflt=0,pos=17)
    review = tables.Int8Col(dflt=0,pos=18)
    cluster_parent = tables.Int64Col(dflt=0,pos=19)


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
    mag_g_wide_limit = tables.Float32Col(dflt=G.HETDEX_CONTINUUM_MAG_LIMIT)
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


class CalibratedSpectra(tables.IsDescription):
    detectid = tables.Int64Col(pos=0)  # unique HETDEX detection ID 1e9+
    #wavelength = tables.Float32Col(shape=(1036,),pos=1) #skip it
    flux = tables.Float16Col(shape=(1036,),pos=2 )  #from Float32Col
    flux_err = tables.Float16Col(shape=(1036,),pos=3)  #from Float32Col
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
    image_depth_mag = tables.Float32Col(dflt=99.9)
    pixel_scale = tables.Float32Col(dflt=UNSET_FLOAT)
    radius = tables.Float32Col(dflt=UNSET_FLOAT) #in arcsec , #was aperture_radius
    mag = tables.Float32Col(dflt=UNSET_FLOAT) #was aperture_mag
    mag_err = tables.Float32Col(dflt=UNSET_FLOAT) #was  aperture_mag_err

    mag_dered = tables.Float32Col(dflt=UNSET_FLOAT) # #added 0.9.2

    aperture_area_pix = tables.Float32Col(dflt=UNSET_FLOAT) #pixels
    sky_area_pix = tables.Float32Col(dflt=UNSET_FLOAT) #pixels
    eqw_rest_lya = tables.Float32Col(dflt=UNSET_FLOAT) #was  aperture_eqw_rest_lya
    eqw_rest_lya_err = tables.Float32Col(dflt=UNSET_FLOAT) #was  aperture_eqw_rest_lya_err
    plae = tables.Float32Col(dflt=UNSET_NAN) #was  aperture_plae
    plae_max = tables.Float32Col(dflt=UNSET_NAN) #was  aperture_plae_max
    plae_min = tables.Float32Col(dflt=UNSET_NAN) #was  aperture_plae_min
    aperture_cts = tables.Float32Col(dflt=UNSET_FLOAT) #was aperture_counts
    sky_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    sky_average = tables.Float32Col(dflt=UNSET_FLOAT)

class ElixerApertures(tables.IsDescription): #perhaps keep just the selected ones??
    #one entry per aperture photometry collected
    detectid = tables.Int64Col(pos=0)
    ra = tables.Float32Col(pos=1,dflt=UNSET_FLOAT) #decimal degrees of center
    dec = tables.Float32Col(pos=2,dflt=UNSET_FLOAT)
    catalog_name = tables.StringCol(itemsize=16)
    filter_name = tables.StringCol(itemsize=16)
    image_depth_mag = tables.Float32Col(dflt=99.9)
    pixel_scale = tables.Float32Col(dflt=UNSET_FLOAT) #arcsec/pixel
    selected = tables.BoolCol(dflt=False) #if True this is the object used for the aperture PLAE/OII, etc (see above table)
    radius = tables.Float32Col(dflt=UNSET_FLOAT) #major axis (diameter) 'a' in arcsec
    mag = tables.Float32Col(dflt=UNSET_FLOAT)
    mag_err = tables.Float32Col(dflt=UNSET_FLOAT)

    mag_dered = tables.Float32Col(dflt=UNSET_FLOAT)  #added 0.9.2

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
    image_depth_mag = tables.Float32Col(dflt=99.9)
    pixel_scale = tables.Float32Col(dflt=UNSET_FLOAT) #arcsec/pixel
    selected = tables.BoolCol(dflt=False) #if True this is the object used for the aperture PLAE/OII, etc (see above table)
    major = tables.Float32Col(dflt=UNSET_FLOAT) #major axis (diameter) 'a' in arcsec
    minor = tables.Float32Col(dflt=UNSET_FLOAT) #'b'
    theta = tables.Float32Col(dflt=0.0) #radians counter-clockwise from x-axis
    mag = tables.Float32Col(dflt=UNSET_FLOAT)
    mag_err = tables.Float32Col(dflt=UNSET_FLOAT)

    mag_dered = tables.Float32Col(dflt=UNSET_FLOAT)  #added 0.9.2


    background_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    background_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_cts = tables.Float32Col(dflt=UNSET_FLOAT)
    flux_err = tables.Float32Col(dflt=UNSET_FLOAT)
    flags = tables.Int32Col(dflt=0)
    dist_curve = tables.Float32Col(dflt=UNSET_FLOAT)
    dist_baryctr = tables.Float32Col(dflt=UNSET_FLOAT)
    image_flags = tables.Int64Col(dflt=0) #separate from the aperture flags, these are ties to the image reduction pipeline

    fixed_aper_radius = tables.Float32Col(dflt=UNSET_FLOAT)
    fixed_aper_mag = tables.Float32Col(dflt=UNSET_FLOAT)
    fixed_aper_mag_err = tables.Float32Col(dflt=UNSET_FLOAT)

    fixed_aper_mag_dered = tables.Float32Col(dflt=UNSET_FLOAT)  #added 0.9.2


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
    match_num = tables.Int32Col(dflt=-1)
    separation = tables.Float32Col(dflt=UNSET_FLOAT) #in arcsec
    prob_match = tables.Float32Col(dflt=UNSET_FLOAT) #in arcsec
    specz = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_specz
    photz = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_photz
    flux = tables.Float32Col(dflt=UNSET_FLOAT) #was  cat_flux
    flux_err = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_flux_err
    mag = tables.Float32Col(dflt=UNSET_FLOAT) #was  cat_mag
    mag_err = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_mag_err
    eqw_rest_lya = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_eqw_rest_lya
    eqw_rest_lya_err = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_eqw_rest_lya_err
    plae = tables.Float32Col(dflt=UNSET_FLOAT) #was  cat_plae
    plae_max = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_plae_max
    plae_min = tables.Float32Col(dflt=UNSET_FLOAT) #was cat_plae_min

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





class VIRUSFiber(tables.IsDescription):
    """
    cloned from HETDEX_API create_shot_hdf5.py

    remove redundant columns (for this purpose ... single shot only)
    reduce the unneccesary precision to save space

    """
    #obsind = tables.Int32Col() #usuall only 3 digits, but I suppose technically could be up to 9, so leave as is?
                                #BUT this is redundant with the shot table since this is just one shot
                                #also redundant with the multiframe string
    multiframe = tables.StringCol((20), pos=0)
    fiber_id = tables.StringCol((38), pos=4)
    # fibidx = tables.Int32Col() #if this is for a single shot, the fibidx is unneccesary
    #                            #as is the finer index table
    fibnum = tables.Int16Col() # only runs 0(1) to 112 or 0(1) to 448
    ifux = tables.Float32Col()
    ifuy = tables.Float32Col()
    fpx = tables.Float32Col()
    fpy = tables.Float32Col()
    ra = tables.Float32Col(pos=1) #may still want this precision, float16 really only gives 4 decimals here
    dec = tables.Float32Col(pos=2)

    calfib = tables.Float16Col((1036,))
    calfibe = tables.Float16Col((1036,))
    calfib_ffsky = tables.Float16Col((1036,)) #could consider dropping this tone too

    ifuslot = tables.StringCol(3)
    ifuid = tables.StringCol(3)
    specid = tables.StringCol(3)
    contid = tables.StringCol(8)
    amp = tables.StringCol(2)
    expnum = tables.Int32Col()

    #consider removing these to save longterm storage
    #if needed, would re-run
    # spectrum = tables.Float16ColCol((1032,))
    # wavelength = tables.Float16ColCol((1032,))
    # fiber_to_fiber = tables.Float16ColCol((1032,))
    #
    # chi2 = tables.Float16ColCol((1032,))
    # rms = tables.Float16ColCol((1032,))
    # calfib_counts = tables.Float16ColCol((1036,))
    # calfibe_counts = tables.Float16ColCol((1036,))
    #
    # trace = tables.Float16ColCol((1032,))
    # sky_subtracted = tables.Float16ColCol((1032,))
    # sky_spectrum = tables.Float16ColCol((1032,))
    # error1D = tables.Float16ColCol((1032,))



def flush_all(fileh,reindex=True):
    # iterate over all tables and issue flush

    if fileh is not None:
        #elixer h5 carry overs
        vtb = fileh.root.Version
        dtb = fileh.root.Detections
        ltb = fileh.root.SpectraLines
        stb = fileh.root.CalibratedSpectra
        atb = fileh.root.Aperture
        ctb = fileh.root.CatalogMatch
        xtb = fileh.root.ElixerApertures
        etb = fileh.root.ExtractedObjects


        vtb.flush()
        dtb.flush()
        ltb.flush()
        stb.flush()
        atb.flush()
        ctb.flush()
        xtb.flush()
        etb.flush()

        #shot h5 carry overs
        shtb = fileh.root.Shot
        ftb = fileh.root.Data.Fibers #in shot h5, this is root.Data.Fibers

        shtb.flush()
        ftb.flush()


        # try: #this table is not always create, so may not exist
        #     ntb = fileh.root.NeighborSpectra
        #     ntb.flush()
        # except:
        #     ntb = None
        #
        # try: #this table is not always created (only with --LyC or --deblend)
        #     dstb = fileh.root.DeblendedSpectra
        #     dstb.flush()
        # except:
        #     dstb = None



        if not reindex:
            return #we're done

        #remove (old) index if exists
        #vtb does not have or need an index
        try:
            dtb.cols.detectid.remove_index()
        except:
            log.debug("Failed to remove detectid index on detections table",exc_info=True)

        try:
            dtb.cols.ra.remove_index()
            dtb.cols.dec.remove_index()
        except:
            log.debug("Failed to remove ra, dec index on detections table",exc_info=True)

        try:
            ltb.cols.detectid.remove_index()
            stb.cols.detectid.remove_index()
            atb.cols.detectid.remove_index()
            ctb.cols.detectid.remove_index()
        except:
            log.debug("Failed to remove detectid index on multiple tables",exc_info=True)

        try:
            ftb.cols.multiframe.remove_index()
            ftb.cols.ra.remove_index()
            ftb.cols.dec.remove_index()
        except:
            log.debug("Failed to remove index on fibers table",exc_info=True)


        dtb.flush()
        ltb.flush()
        stb.flush()
        atb.flush()
        ctb.flush()

        ftb.flush()

        #create (new) index
        # vtb does not have or need an index
        try:
            dtb.cols.detectid.create_csindex()
        except:
            log.debug("Index fail on detections table: detectid",exc_info=True)

        try:
            dtb.cols.ra.create_csindex()
            dtb.cols.dec.create_csindex()
        except:
            log.debug("Index fail on detections table: ra and/or dec",exc_info=True)


        try:
            ltb.cols.detectid.create_csindex()
        except:
            log.debug("Index fail on lines table",exc_info=True)

        try:
            stb.cols.detectid.create_csindex()
        except:
            log.debug("Index fail on spectra table",exc_info=True)

        try:
            atb.cols.detectid.create_csindex()
        except:
            log.debug("Index fail on apertures table",exc_info=True)

        try:
            ctb.cols.detectid.create_csindex()
        except:
            log.debug("Index fail on catalog match table",exc_info=True)


        #vtb.flush() # no need to re-flush vtb
        dtb.flush()
        ltb.flush()
        stb.flush()
        atb.flush()
        ctb.flush()

        try:
            etb.cols.detectid.remove_index()
            etb.cols.detectid.create_csindex()
            etb.flush()
        except:
            log.debug("Index fail on extracted objects table")

        try:
            xtb.cols.detectid.remove_index()
            xtb.cols.detectid.create_csindex()
            xtb.flush()
        except:
            log.debug("Index fail on elixer apertures table")

        try:
            ftb.cols.multiframe.create_csindex()
            ftb.cols.ra.create_csindex()
            ftb.cols.dec.create_csindex()
        except:
            log.debug("Index fail on Fibers table",exc_info=True)


        ftb.flush()

        #
        # try:
        #     if ntb is not None:
        #         ntb.cols.detectid.remove_index()
        #         ntb.cols.detectid.create_csindex()
        #         ntb.flush()
        # except:
        #     log.debug("Index fail on neighbor spectra table")
        #
        # try:
        #     if dstb is not None:
        #         dstb.cols.detectid.remove_index()
        #         dstb.cols.detectid.create_csindex()
        #         dstb.flush()
        # except:
        #     log.debug("Index fail on deblended spectra table")
        #
        # try:
        #     if vote_tb is not None:
        #         vote_tb.cols.detectid.remove_index()
        #         vote_tb.cols.detectid.create_csindex()
        #         vote_tb.flush()
        # except:
        #     log.debug("Index fail on vote table")

    return


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


    try:
        shot_h5 = None
        elixer_h5 = None

        #check that the one or two input files exist and can be opened
        if shot_fn is not None:
            try:
                shot_h5 = tables.open_file(shot_fn,mode='r')
            except:
                log.critical(f"Failed to open input shot.h5 file: {shot_fn}",traceback.format_exc())
                return

        if shot_h5 is None: #this MUST always be provided
            log.critical(f"Fatal. Cannot open shot.h5 file: {shot_fn}")
            return

        if elixer_fn is not None:
            try:
                elixer_h5 = tables.open_file(elixer_fn,mode='r')
            except:
                log.critical(f"Fatal. Failed to open input elixer.h5 file: {elixer_fn}",traceback.format_exc())
                return

            if elixer_h5 is None:
                log.critical(f"Fatal. Failed to open input elixer.h5 file: {elixer_fn}", traceback.format_exc())
                return

            #check supported version
            elixer_h5_version = elixer_h5.root.Version.read(field='version')[0]
            if elixer_h5_version not in SUPPORTED_ELIXER_H5_VERSIONS:
                log.critical(f"Fatal. Input elixer h5 version {elixer_h5_version} not supported.")
                return


        #todo: ?? should we name differently if elixer is present vs not?
        # say "r"+shot_fn if just the shot
        # and "re" + shot_fn if shot and elixer
        # (or "s" and "se") ??


        outfn = "ssr" + os.path.basename(shot_fn)

        log.debug("Creating new SingleShot Reduction HDF5 catalog (%s)" % (outfn))

        fileh = tables.open_file(outfn, 'w', 'SingleShot Reduction Catalog')

        vtb = fileh.create_table(fileh.root, 'Version', Version,
                                 'Version Table')

        row = vtb.row
        row['version'] = __version__
        row['version_pytables'] = tables.__version__
        row.append()
        vtb.flush()

        fileh.create_table(fileh.root, 'Shot', VIRUSShot,
                           'Shot Summary Table')

        #put under "Data" for some compatibility with pathing in HETDEX_API
        group_data = fileh.create_group(fileh.root, "Data", "VIRUS Fiber Data")
        fileh.create_table(group_data, 'Fibers', VIRUSFiber, 'Fiber Summary Table')

        #iterate over shot's info (shot table and fiber table) and populate

        #just one row for Shot table, so just do a direct copy
        fileh.root.Shot.append(shot_h5.root.Shot.read())
        fileh.root.Shot.flush()

        #many for fibers, so iterate
        for row in shot_h5.root.Data.Fibers.read():
            new_row = fileh.root.Data.Fibers.row
            #go over some columns individual since want to change some types
            #directy copy columns:
            for col in ['multiframe','fiber_id','ifux','ifuy','fpx','fpy','ra','dec',
                        'ifuslot','ifuid','specid','contid','amp']:
                new_row[col] = row[col]

            #special treatment, mostly float32 to float16
            #expnum (int32 to int16)
            #fibnum (int32 to int16)
            #calfib (float32 array to float16 array) #Float16Col((1036,))
            #calfibe (float32 array to float16 array) #Float16Col((1036,))
            #calfib_ffsky (float32 array to float16 array) #Float16Col((1036,))

            new_row['expnum'] = row['expnum'].astype(np.int16) #Maybe we'd have more than 15 exposures (prob not, though)
            new_row['fibnum'] = row['fibnum'].astype(np.int16)
            new_row['calfib'] = row['calfib'].astype(np.float16)
            new_row['calfibe'] = row['calfibe'].astype(np.float16)
            new_row['calfib_ffsky'] = row['calfib_ffsky'].astype(np.float16)

            new_row.append()

        fileh.root.Data.Fibers.flush()
        #set the index
        try:
            fileh.root.Data.Fibers.cols.multiframe.create_csindex()
            fileh.root.Data.Fibers.cols.ra.create_csindex()
            fileh.root.Data.Fibers.cols.dec.create_csindex()
        except:
            log.debug("Index fail on fibers table", traceback.format_exc())

        fileh.root.Data.Fibers.flush()
        #done with the shot.h5 file now
        shot_h5.close()



        if elixer_h5 is not None:

            fileh.create_table(fileh.root, 'Detections', Detections,
                               'Detection Summary Table')

            fileh.create_table(fileh.root, 'SpectraLines', SpectraLines,
                               'Identified SpectraLines Table')

            fileh.create_table(fileh.root, 'CalibratedSpectra', CalibratedSpectra,
                               'PSF Weighted Spectra Table')

            fileh.create_table(fileh.root, 'Aperture', Aperture,
                               'Aperture Photometry Table') # mostly a g and r aperture, sometimes more

            fileh.create_table(fileh.root, 'CatalogMatch', CatalogMatch,
                               'Archival Catalog Matched Objected Table')

            fileh.create_table(fileh.root, 'ExtractedObjects', ExtractedObjects,
                               'SourceExtractor Objects Table')  # multiple filters, many objects

            fileh.create_table(fileh.root, 'ElixerApertures', ElixerApertures,
                               'Circular Apertures Table')  # mostly a g and r aperture, sometimes more

            ###########################################################################
            #Detections all direct copy but not all columns
            ###########################################################################
            for row in elixer_h5.root.Detections.read():
                new_row = fileh.root.Detections.row

                # go over some columns individual since want to change some types
                # directy copy columns:

                skip_cols = np.array(['shotid','obsid','seeing_fwhm','response','fieldname'])
                keep_cols = np.setdiff1d(elixer_h5.root.Detections.colnames,skip_cols)

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
                log.debug("Index fail on Detections table", traceback.format_exc())

            fileh.root.Detections.flush()

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
                log.debug("Index fail on SpectraLines table", traceback.format_exc())

            fileh.root.SpectraLines.flush()

            #############################################################################
            # CalibratedSpectra
            #############################################################################
            for row in elixer_h5.root.CalibratedSpectra.read():
                new_row = fileh.root.CalibratedSpectra.row

                new_row['detectid'] = row['detectid']
                # skip wavelength entirely and set the float32 to float16
                new_row['flux'] = row['flux'].astype(np.float16)
                new_row['flux_err'] = row['flux_err'].astype(np.float16)
                new_row['dust_corr'] = row['dust_corr'].astype(np.float16)
                new_row['aperture_radius'] = row['aperture_radius'].astype(np.float16)
                new_row['sky_background'] = row['sky_background'].astype(np.int8)
                new_row['num_fibers'] = row['num_fibers'].astype(np.int16)

                new_row.append()

            fileh.root.CalibratedSpectra.flush()
            #set the index
            try:
                fileh.root.CalibratedSpectra.cols.detectid.create_csindex()
            except:
                log.debug("Index fail on CalibratedSpectra table", traceback.format_exc())

            fileh.root.CalibratedSpectra.flush()



            # done with the elixer.h5 file now
            elixer_h5.close()


    except:
        print(f"Exception building new h5 file. {traceback.format_exc()}")



    #!! todo for ELiXerApertures, only keep the "selected==True" rows




    fileh.close()
#end build_ssr_shot_h5



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

if len(args) > 0:
    print(f"Unknown remainting args: {args}")

build_ssr_shot_h5(shot_h5_path, elixer_fn=elixer_h5_path)