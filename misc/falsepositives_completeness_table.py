"""
see also dustin laptop notebook: dets_cross_match-working


this has to be in the directory with the base tables

"""


from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
#import astropy.io.fits as fits
import tables
import glob
import os
import socket
#import copy
from elixer import utilities as U
from elixer import global_config as G
#from hetdex_api.elixer_widget_cls import ElixerWidget
#from hetdex_tools.get_spec import get_spectra as get_spectra
#from hetdex_api.shot import get_image2D_amp
#from hetdex_api.amp_widget import AmpWidget
#from astrowidgets import ImageWidget
#flux limit includes
from hetdex_api.flux_limits import flim_models
from hetdex_api.shot import Fibers
import traceback
from tqdm import tqdm
import numpy as np
import time
from datetime import datetime
#import astropy.stats.biweight as biweight

import sys
SHOW_TQDM = False

if "tacc" in socket.gethostname():
    LOCAL = False
    SHOW_TQDM = False
else:
    LOCAL = True

args = list(sys.argv) #python3 map is no longer a list, so need to cast here
del args[0] #args.pop(0) #remove THIS file
args = [x.replace("--","-") for x in args]

#then single parameter is the group number

group_num = int(args[0])
dtnow = datetime.now()
dtnowstr = "%s:%s:%s.%s : " % (str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
                            str(dtnow.microsecond).zfill(6))
print(f"{dtnowstr} Working on group {str(group_num).zfill(2)}")
u_groups = [group_num]

if not LOCAL and not os.path.exists(f"group_{str(group_num).zfill(2)}"):
    print("No such group. Exiting...")
    exit(-1)

#loop over the shots in a group part


########################################
# need a dictionary of shots:
#######################################
shot_fiber_dict = {}

shot_fns = glob.glob(f"/scratch/03261/polonius/parallel/line_det_test/group_{str(group_num).zfill(2)}/ssr_*.h5")

#
# if LOCAL:
#     #shot_fns = glob.glob("/run/media/dustin/easystore/laurel_repeats/ssr_*.h5")
#     shot_fns = glob.glob(f"/scratch/03261/polonius/parallel/line_det_test/group_{str(group_num).zfill(2)}/ssr_*.h5")
#     #now, get the base name and replace the directory
#
#     shot_fns = [
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250821v019.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250824v022.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250825v021.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250826v017.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250827v009.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250829v019.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250922v019.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20250923v023.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20251017v008.h5",
#         "/run/media/dustin/easystore/laurel_repeats/ssr_20251022v015.h5",
#
#     ]
# else:
#     shot_fns = glob.glob(f"/scratch/03261/polonius/parallel/line_det_test/group_{str(group_num).zfill(2)}/ssr_*.h5")

if len(shot_fns) < 1:
    print("FATAL. Could not find shots.")
    exit(-1)

for fn in shot_fns:
    if LOCAL:
        fn = f"/run/media/dustin/easystore/laurel_repeats/" + os.path.basename(fn)

    print(f"Loading {fn} ...")
    shot = int(os.path.basename(fn).split("ssr_")[1].split(".")[0].replace("v", ""))
    d = {}
    d['Fibers'] = Fibers(shot=shot, shot_h5=fn)

    ssr_h5 = tables.open_file(fn)
    shot_row = ssr_h5.root.Shot.read()
    d['seeing'] = shot_row['fwhm_virus'][0]
    d['tput'] = shot_row['response_4540'][0]
    ssr_h5.close()

    shot_fiber_dict[shot] = d

#####################################
#Load the base tables we need
#####################################

T_mc_dets = Table.read('laurel_mc_dets.fits',format="fits")
T_elixer_dets = Table.read('laurel_elixer_dets.fits',format="fits")
T_fiber_coords = Table.read('laurel_fiber_coords.fits',format="fits")
T_elixer_sep = Table.read('laurel_elixer_sep.fits',format="fits")

if not "shotid" in T_elixer_dets.indices:
    T_elixer_dets.add_index("shotid")
if not "group" in T_elixer_dets.indices:
    T_elixer_dets.add_index("group")

if not "shotid" in T_elixer_sep.indices:
    T_elixer_sep.add_index("shotid")
if not "group" in T_elixer_sep.indices:
    T_elixer_sep.add_index("group")

if not "shotid" in T_fiber_coords.indices:
    T_fiber_coords.add_index("shotid")
if not "group" in T_fiber_coords.indices:
    T_fiber_coords.add_index("group")

if not "shotid" in T_mc_dets.indices:
    T_mc_dets.add_index("shotid")
if not "group" in T_mc_dets.indices:
    T_mc_dets.add_index("group")


fixed_sn_cut = 4.8
PSFtab = Table.read("psf_fiber_convolutions_radius_3.5.txt", delimiter="\t", format='ascii')
PSFcols = np.array(list(PSFtab.columns[1:])).astype(float)


# the noise used in the creation, test, calibration of the f50 stuff is the psf weighted errror
# spectrum_error = np.sqrt(np.sum(error ** 2 * mask * weights, axis=0) / np.sum(mask * weights ** 2, axis=0)
# ... I do think "error" in the func call above is the calfibe ... so it is the normalized, quadrature sum with (PSF) weights
# .... what you get back is an array of "errors" at each wavelength bin in the extraction ....

# we don't really have that at this stage since there is no detection ... we *could* do an extraction AT that position,
# but that is not really quite fair and it is very slow when done out over eveything ...
# we CANNOT use the error from the source detection and cannot use its weights or mask ... we need it from the target bid shots
# ...

# so

# todo LETS DO THE FOLLOWING
# make an array of weights assuming a 3.0 (or 3.5?) arcsec extraction for seeings from 1.2 to 3.0
# then for each of the target (bid) shots that have no detection, grab their calfibe (like we are now), and apply the weights
# (ignore the mask? that is ... assume they are all good? or could get the mask array too and apply it)
# choose the nearest PSF seeing (1.2 to 3.0) for the weights and use the "spectrum_error" compuation above to compute the nosise to
# use for the f50 stuff

# there should be a PSF notebook somewhere I have already done that can be used for the weights ... or just extract the function and run
# for PSF seeing 1.2 to 3.0

def get_psf_conversion(wavelength, target_psf):
    """
    a previous generated table that approximates the conversion from one wavelength bin to a whole PSF
    (e.g. use the returned float to multiply into the value of that wavelenth bin (the flux or error)
    and you get the weighted total for the PSF as IF all the fibers inside the PSF had that same value at their same wavelength)
    see notebook "psf_convolution_keys" which builds the psf_fiber_convolutions_radius_3.5.txt file
    """

    ix, *_ = U.getnearpos(PSFcols, target_psf)
    col = f"{PSFcols[ix]:0.2f}"
    ix, *_ = U.getnearpos(G.CALFIB_WAVEGRID, wavelength)
    row = int(G.CALFIB_WAVEGRID[ix])

    return float(PSFtab[col][PSFtab['wave'] == row])


def get_noise(shot_fibers_obj, ra, dec, radius, source_sigma, source_wavelength, target_psf):
    try:
        num_sigma = 2.0
        num_fibers = -1
        # _,left,_ = U.getnearpos(G.CALFIB_WAVEGRID,source_wavelength-num_sigma*source_sigma)
        # _,_,right = U.getnearpos(G.CALFIB_WAVEGRID,source_wavelength+num_sigma*source_sigma)
        # if right is None:
        #     right = len(G.CALFIB_WAVEGRID)
        # if left is None:
        #     left = 0

        # just use the closest, so don't have to check that we are off one end or the other
        left, *_ = U.getnearpos(G.CALFIB_WAVEGRID, source_wavelength - num_sigma * source_sigma)
        right, *_ = U.getnearpos(G.CALFIB_WAVEGRID, source_wavelength + num_sigma * source_sigma)

        ff = shot_fibers_obj.query_region(SkyCoord(ra=ra * u.deg, dec=dec * u.deg), radius=radius)
        num_fibers = len(ff)

        # near the edges, the error can be 0.0 and should NOT be used
        ff['calfibe'][ff['calfibe'] == 0.0] = np.nan

        avg_calfibe = np.nanmedian(ff['calfibe'][:, left:right + 1])
        psf_col = get_psf_conversion(source_wavelength, target_psf)
        # last remaining question is whether I need the 2x for flux from fluxd
        # right now, I think so, since calfibe, as read in, is fluxd
        # PLUS at Seeing FWHM 1.5, this agrees pretty well with Karl's plot of completness in the paper
        avg = avg_calfibe * psf_col * 2.0

        # the average of all the fibers withing radius (usually 3") for the wavelength bins within +/- (2.5?) sigma of line center
        # this is in per AA, so then multiply by 2 * 2.5 * sigma to get the AA coverage

        # should this weight better? since it would be a PSF extraction
        # instead of 2 * num_sigma * source_sigma, use the Gaussian area? (since these are flux units, the Area is the flux)
        # so could replace the 2*num_sigma*source_sigma with * lineflux) -- NO !!! ... while the flux is an "area" it is not
        # a weighting solution (think about the same linewidth but twice as much lineflux ... if using the lineflux as a weight
        # you'd get a much worse completness since the "noise" would be twice as high ... the exact opposite of what should be
        #  (twice the flux would be a larger completeness)

        #   +/- 1 sigma  * 2AA/bin * 1/0.68 = 2.94
        # or +/- 2 sigma  * 2AA/bin * 1/0.95 = 4.21
        # avg = 4.21 * source_sigma * np.nanmedian(ff['calfibe'][:,left:right+1])
        # naive form
        # avg = 2 * num_sigma * source_sigma * np.nanmedian(ff['calfibe'][:,left:right+1])

        # avg = (right - left + 1) * np.nanmedian(ff['calfibe'][:,left:right+1])
        # print("todo: check is this error for flux density or flux ?")
        # these are flux density units ... so, per AA ... so depending on how we count, we need to * 2
        # avg = 2.0 * np.nanmax(ff['calfibe'][:,left:right+1]) #per AA, so need * 2
    except:
        print(traceback.format_exc())
        avg = None
        num_fibers = -1

    return avg, num_fibers


def get_f50(shot_fibers_obj, ra, dec, radius, source_sigma, source_wavelength, target_psf):
    try:
        noise, num_fibers = get_noise(shot_fibers_obj, ra, dec, radius, source_sigma, source_wavelength, target_psf)
        f50 = flim_models.f50_from_noise_v5_simple(noise=noise, wave=source_wavelength, sncut=fixed_sn_cut,
                                                   linewidth=source_sigma)
        f50_div2 = f50  # flim_models.f50_from_noise_v5_simple(noise=noise/2.0, wave=source_wavelength, sncut=fixed_sn_cut,linewidth=source_sigma)
    except:
        print(traceback.format_exc())
        f50 = None
        f50_div2 = None
        num_fibers = None

    return f50, f50_div2, num_fibers


# def get_f50(shot_fibers_obj,ra,dec,radius,source_sigma,source_wavelength,source_fluxe=0.0,):
#     try:
#         if source_fluxe==0.0:
#             noise = get_noise(shot_fibers_obj,ra,dec,radius,source_sigma,source_wavelength) #binning
#         else:
#             noise = source_fluxe ** 2
#         f50 = flim_models.f50_from_noise_v5_simple(noise=noise, wave=source_wavelength, sncut=fixed_sn_cut,linewidth=source_sigma)
#     except:
#         print(traceback.format_exc())
#         f50 = None

#     return f50

def get_limits(shot_fibers_obj, ra, dec, radius, source_sigma, source_wavelength, source_flux, target_psf):
    """
    return the completeness and f50 limit
    source_flux in e-17 scale
    source_sigma in AA
    source_wavelength in AA
    """

    if isinstance(shot_fibers_obj, dict):
        shot_fibers_obj = shot_fibers_obj['Fibers']
    f50, f50_div2, num_fibers = get_f50(shot_fibers_obj, ra, dec, radius, source_sigma, source_wavelength, target_psf)
    # print("f50",f50)
    cmplt = flim_models.completeness_v5(flux=source_flux, f50=f50, wave=source_wavelength, sncut=fixed_sn_cut)
    cmplt_div2 = cmplt  # flim_models.completeness_v5(flux=source_flux, f50=f50_div2, wave=source_wavelength, sncut=fixed_sn_cut)
    return f50, cmplt, f50_div2, cmplt_div2, num_fibers





################################
#only send in one group
###############################

if True: #response.lower() == 'y':

    skip_stats_on_match_found = False  # if true, do not execute the section that collects the f50, compeleteness when a match WAS FOUND

    source_shotid = []
    source_group = []
    source_detid = []
    miss_shotid = []
    # source_flux = []
    # miss_seeing = []
    # miss_avg_calfibe = []
    miss_f50 = []
    miss_completeness = []

    miss_f50_2 = []
    miss_completeness_2 = []
    miss_num_fibers = []

    # miss_apcor = [] #unless we re-extract, don't get this one
    miss_seeing = []
    miss_tput = []
    miss_found = []  # was this actually found (not really a miss)
    miss_self = []  # T/F is this the same as the source shot
    miss_sn = []

    # repeat_found = []
    source_seeing = []
    source_tput = []
    source_wave = []
    source_line_flux = []
    source_line_sigma = []

    source_sn = []  # new

    source_f50 = []
    source_completeness = []

    source_f50_2 = []
    source_completeness_2 = []
    source_num_fibers = []

    halt = False

    # out_dets = []
    # out_shots = [] #for future calls
    # out_sn = []
    # out_matches = []
    # out_miss_ct = []
    # out_has_multiples = [] #e.g. in a single shot ... put in the list of shotids

    # to be considered a match
    base_max_sep = 2.0  # arcsec  ... this actually needs to be flexible, esp for resolved objects !!!
    base_max_resolved_major = 15.0 #arcsec
    max_delta_wave = 6.0  # AA   #this is the usual, but I worry it gets a bit to lax for small lin ewidths

    # no-match found check on fibers
    aperture = 3.5  # arcsec ... this is the aperture using the target coords on the bid_shot, assume some
    # astrometric error between the target and bid shots, so make larger than the normal 3.0 - 3.5"
    # will not be quite right for extended objects, but those are handled a little differently
    min_fiber_ct = 9  # how many fibers should be in the aperture (e.g. this keeps us away from the edge)
    # and keeps us from counting as misses the shots that don't have good matching fiber coverage

    # reset?
    T_elixer_dets['num_matches'] = -1
    T_elixer_dets[
        'num_valid_misses'] = -1  # in how many group shots was it not found, where there are sufficient fibers

    already_done_shots = [] #empty list, but read in those already done if table found on disk

    #u_groups = np.unique(T_elixer_dets['group'])

    # #need a dictionary of shots:
    # print("Loading fiber dicts for shots ...")
    # shot_fiber_dict = {}
    # shot_fns = glob.glob("/scratch/03261/polonius/parallel/line_det_test/group_??/ssr_*.h5")
    # for fn in tqdm(shot_fns):
    #     shot = int(os.path.basename(fn).split("ssr_")[1].split(".")[0].replace("v",""))
    #     shot_fiber_dict[shot] = Fibers(shot=shot,shot_h5=fn)


    # print("just running 1 group, 1 shot")
    for group in u_groups:  # [0:1]):
        dtnow = datetime.now()
        dtnowstr = "%s:%s:%s.%s : " % (str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
                                    str(dtnow.microsecond).zfill(6))
        print(f"{dtnowstr} Starting group {group}")
        # if halt:
        #     break
        # get the group sub-selected tables
        td = T_elixer_dets[T_elixer_dets['group'] == group]
        tf = T_fiber_coords[T_fiber_coords['group'] == group]
        Te = T_elixer_sep[T_elixer_sep['group'] == group]
        tmX = T_mc_dets[T_mc_dets['group'] == group]

        ra_for_dec_adjust = np.cos(np.deg2rad(td['dec'][0]))  # degrees, adjusted for declination

        # need to compute for EACH detection as the size (e.g. if it is resovled) matters
        # max_dec = max_sep / 3600.0 #in degrees
        # max_ra = max_dec / ra_for_dec_adjust
        # #all groups are similar in dec and all detects are similar (with 15' or so), so this is close enough
        # max_miss_dec = (aperture + 1.5) / 3600.0 #go out larger than the aperture
        # max_miss_ra = max_miss_dec / np.cos(np.deg2rad(td['dec'][0]))

        # iterate over each shot
        u_shots = np.unique(td['shotid'])

        # print(f"group {group}, {len(u_shots)} shots)")

        #if the table is already there, read it in and resume from the last write
        if os.path.exists(f"laurel_elixer_mc_completeness_{str(group_num).zfill(2)}.fits"):
            T = Table.read(f"laurel_elixer_mc_completeness_{str(group_num).zfill(2)}.fits", format="fits")

            source_group = list(T['group'])
            source_shotid = list(T['src_shotid'])
            source_detid = list(T['detectid'])
            source_wave = list(T['wave'])
            source_line_flux = list(T['line_flux'])
            source_line_sigma = list(T['line_sigma'])
            source_num_fibers = list(T['src_nfib'])
            source_sn = list(T['sn'])
            source_f50 = list(T['src_f50'])
            source_completeness = list(T['src_cmplt'])
            miss_shotid = list(T['bid_shotid'])
            miss_self = list(T['self'])
            miss_found = list(T['found'])
            miss_sn = list(T['bid_sn'])
            miss_seeing = list(T['bid_seeing'])
            miss_tput = list(T['bid_tput'])
            miss_num_fibers = list(T['bid_nfib'])
            miss_f50 = list(T['bid_f50'])
            miss_completeness = list(T['bid_cmplt'])

            already_done_shots = np.unique('source_shotid')



        # for i, zp in enumerate(tqdm(td)):
        shot_counter = 0
        for shot in tqdm(u_shots,disable=not SHOW_TQDM):  # [0:1]):
            # if halt:
            #     break

            shot_counter += 1
            dtnow = datetime.now()
            dtnowstr = "%s:%s:%s.%s : " % (str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
                                str(dtnow.microsecond).zfill(6))
            if shot in already_done_shots:
                print(f"{dtnowstr} Already done group {group} shot {shot_counter} / {len(u_shots)} : {shot} with {len(elixer_dets)} dets ...")
                continue
            else:
                print(f"{dtnowstr} Starting group {group} shot {shot_counter} / {len(u_shots)} : {shot} with {len(elixer_dets)} dets ...")


            elixer_dets = td[td['shotid'] == shot]['detectid']
            # print(f"#dets {len(elixer_dets)}")
            # T =zp[1]

            for det in tqdm(elixer_dets,disable=not SHOW_TQDM):  # [0:200]):
                # if halt:
                #     break

                # we dont want continuum sources
                if str(det)[9] == '9':
                    continue

                max_sep = base_max_sep  # reset for each detection
                found_but_get_stats_anyway = False
                repeat_hits = 0
                # out_dets.append(det)
                # out_matches.append([])
                # out_miss_ct.append(0)
                # out_has_multiples.append([])

                sel_update = T_elixer_dets['detectid'] == det
                # ri = T_elixer_dets[T_elixer_dets['detectid']==det]
                ri = T_elixer_dets[sel_update]

                # out_sn.append(ri['sn'][0])
                # out_shots.append(ri['shotid'][0])

                # need the size of the detection ... the RA, Dec differences need to be allowed to occur within the
                # SEP ellipse
                try:
                    sel_te = Te['detectid'] == det
                    if np.count_nonzero(sel_te) > 0:
                        te = Table(Te[sel_te])
                        sel_te = te['selected'] == True
                        major = np.nanmedian(te['major'][sel_te])
                        if not np.isnan(major):
                            max_sep = min(base_max_resolved_major,max(max_sep, major))
                            # if major > max_sep:
                            #     #print(f"Changed max_sep for major = {major}")
                            #     max_sep = major
                        else:
                            print(f"major was NaN; keeping max_sep {max_sep}")
                except:
                    print(f"Exception with {det} for SEP")

                # max_dec = max_sep / 3600.0 #in degrees
                # max_ra = max_dec / ra_for_dec_adjust #set earlier
                # all groups are similar in dec and all detects are similar (with 15' or so), so this is close enough
                #max_miss_dec = (max_sep + 1.5) / 3600.0  # go out a little larger than the aperture
                max_miss_dec = max_sep  / 3600.0  #do not add anything extra ... just build it into the base
                max_miss_ra = max_miss_dec / np.cos(np.deg2rad(td['dec'][0]))

                ref_coord = SkyCoord(ra=ri['ra'] * u.deg, dec=ri['dec'] * u.deg)

                # is this a bad detection?? if so skip
                # FUTURE can do better with new reductions and the full fiber flagging is available per wavelength
                # and could just look at the wavelength range

                # technically we *could* have this off for extended objects, but that is rare enough to not worry about
                # if an extended object falls off the IFU where we want to look, then this just does not add to the
                # set of valid shots to check ... could have a small impact to either the confirmation or failure side, but
                # it just won't count toward either
                ff = shot_fiber_dict[shot]['Fibers'].query_region(ref_coord, radius=max(aperture, max_sep / 2.0))
                if np.any(ff['flag'] == 0):  # this is bad in some way
                    print(f"[{det}] bad ff")
                    continue

                # full length of td
                allowed_delta_wave = min(2.5 * float(ri['fwhm_line_aa']) / 2.355, max_delta_wave)
                sel_wave = np.array(abs(td['wavelength_obs'] - ri['wavelength_obs']) <= allowed_delta_wave)

                # td_bid_match = td[sel_wave]
                sel_radec = np.array(abs(td['ra'][sel_wave] - ri['ra']) <= max_miss_ra) * np.array(
                    abs(td['dec'][sel_wave] - ri['dec']) <= max_miss_dec)

                # includes itself (this is a longer calculation ... first check wave)
                # only length of sub-selection (of the sel_wave that are True)
                bid_det_coords = SkyCoord(ra=np.array(td['ra'][sel_wave][sel_radec]) * u.deg,
                                          dec=np.array(td['dec'][sel_wave][sel_radec]) * u.deg)

                sep = ref_coord.separation(bid_det_coords)
                sel_coord = sep.arcsec <= max_sep

                # so now how to combine
                # td_bid_match = td[sel_wave][sel_radec][sel_coord]

                matching_shots = np.unique(td['shotid'][sel_wave][sel_radec][sel_coord])
                ct = len(matching_shots)  # only count a shot once, even if more than
                # one match is found (rare, but does happen)

                if ct > 0:
                    x_match_shot_list = td['shotid'][sel_wave][sel_radec][sel_coord]
                    x_match_sn_list = td['sn'][sel_wave][sel_radec][sel_coord]
                    match_shot_list = []
                    match_sn_list = []
                    for st in matching_shots:
                        sst = x_match_shot_list == st
                        match_shot_list.append(st)
                        match_sn_list.append(np.max(x_match_sn_list[sst]))
                else:
                    match_shot_list = []
                    match_sn_list = []

                # need this from the source for later
                d = shot_fiber_dict[shot]
                src_seeing = d['seeing']
                src_tput = d['tput']
                src_f50, src_completeness, src_f50_2, src_completeness_2, src_num_fibers = get_limits(
                    d['Fibers'], float(ri['ra']), float(ri['dec']),
                    3.5, float(ri['fwhm_line_aa']) / 2.355,
                    float(ri["wavelength_obs"]), float(ri["flux_line"]) * 1e17, d["seeing"])

                if ct < 1:  # there is a problem ... it should have found itself
                    print(f"Problem: {det} did not find itself.")
                    continue  # leave the matches at -1
                else:  # ct is 1 or greater
                    # if (ct > 1):
                    #    print(f"[{det}] found : {ct}")
                    valid_misses = 0
                    mc_hits = 0

                    for bid_shot in u_shots:
                        # don't skip we want to record it too
                        # if bid_shot == shot:
                        #     continue #skip itself

                        sel_tf = tf['shotid'] == bid_shot
                        sel_radec = np.array(abs(tf['ra'][sel_tf] - ri['ra']) <= max_miss_ra) * \
                                    np.array(abs(tf['dec'][sel_tf] - ri['dec']) <= max_miss_dec)

                        bid_fiber_coords = SkyCoord(ra=np.array(tf['ra'][sel_tf][sel_radec]) * u.deg,
                                                    dec=np.array(tf['dec'][sel_tf][sel_radec]) * u.deg)

                        sep = ref_coord.separation(bid_fiber_coords)

                        num_in_range = np.count_nonzero(sep.arcsec < aperture)
                        # collect the completeness info regardless of whether it was found or not OR is THIS shot

                        okay_to_continue = False
                        if bid_shot == shot:  # it is the original source shot
                            okay_to_continue = True
                        elif num_in_range >= min_fiber_ct:
                            # is this a bad (bid) detection?? if so skip
                            # FUTURE can do better with new reductions and the full fiber flagging is available per wavelength
                            # and could just look at the wavelength range
                            ff = shot_fiber_dict[bid_shot]['Fibers'].query_region(ref_coord,
                                                                                  radius=max(aperture, max_sep / 2.0))
                            if not np.any(ff['flag'] == 0):  # this not explicity bad
                                okay_to_continue = True

                        if okay_to_continue:  # assume coverage is sufficient
                            # now ... check those *.mc files ... could be match, just without sufficient SNR

                            if bid_shot not in match_shot_list:  # did not already found a match in ELiXer

                                # want the mc table selection that matches bid_shot
                                # tm = T_mc_dets[T_mc_dets['shotid'] == bid_shot] #shotid is only in any one group
                                tm = tmX[tmX['shotid'] == bid_shot]  # shotid is only in any one group

                                # full length of td
                                sel_mc_wave = np.array(abs(tm['wave'] - ri['wavelength_obs']) <= max_delta_wave)

                                # td_bid_match = td[sel_wave]
                                sel_mc_radec = np.array(
                                    abs(tm['ra'][sel_mc_wave] - ri['ra']) <= max_miss_ra) * np.array(
                                    abs(tm['dec'][sel_mc_wave] - ri['dec']) <= max_miss_dec)

                                # includes itself (this is a longer calculation ... first check wave)
                                # only length of sub-selection (of the sel_wave that are True)
                                bid_det_mc_coords = SkyCoord(ra=np.array(tm['ra'][sel_mc_wave][sel_mc_radec]) * u.deg,
                                                             dec=np.array(tm['dec'][sel_mc_wave][sel_mc_radec]) * u.deg)

                                sep_mc = ref_coord.separation(bid_det_mc_coords)
                                sel_mc_coord = sep_mc.arcsec <= max_sep

                                # matched wave and position
                                if np.count_nonzero(sel_mc_coord) > 0:
                                    maxsn = max(tm['sn'][sel_mc_wave][sel_mc_radec][sel_mc_coord])
                                    match_shot_list.append(bid_shot)
                                    match_sn_list.append(maxsn)

                            # get the fiber spectra, but do not re-extract?

                            # todo: need seeing fwhm, etc.

                            d = shot_fiber_dict[bid_shot]
                            try:
                                # miss_seeing.append(d["seeing"])
                                # def get_limits(shot_fibers_obj,ra,dec,radius,source_sigma,source_wavelength,source_flux,target_psf):
                                f50, completeness, f50_2, completeness_2, num_fibers = get_limits(
                                    d['Fibers'], float(ri['ra']), float(ri['dec']),
                                    3.5, float(ri['fwhm_line_aa']) / 2.355,
                                    float(ri["wavelength_obs"]), float(ri["flux_line"]) * 1e17, d["seeing"])

                                # miss_avg_calfibe.append()
                                # repeat_found.append(0)

                                sst = match_shot_list == bid_shot
                                if np.count_nonzero(sst) == 1:
                                    miss_sn.append(float(np.array(match_sn_list)[sst]))
                                else:
                                    miss_sn.append(0.)

                                miss_f50.append(f50)
                                miss_completeness.append(completeness)

                                miss_f50_2.append(f50_2)
                                miss_completeness_2.append(completeness_2)
                                miss_num_fibers.append(num_fibers)

                                source_shotid.append(shot)
                                source_group.append(group)
                                source_detid.append(det)
                                miss_shotid.append(bid_shot)
                                if bid_shot in match_shot_list or bid_shot == shot:
                                    miss_found.append(True)
                                else:
                                    miss_found.append(False)
                                    valid_misses += 1

                                # miss_apcor.append()
                                miss_seeing.append(d["seeing"])
                                miss_tput.append(d["tput"])

                                source_seeing.append(src_seeing)
                                source_tput.append(src_tput)
                                source_wave.append(float(ri["wavelength_obs"]))
                                source_line_flux.append(float(ri["flux_line"]) * 1e17)
                                source_line_sigma.append(float(ri['fwhm_line_aa']) / 2.355)

                                source_f50.append(src_f50)
                                source_completeness.append(src_completeness)

                                source_f50_2.append(src_f50_2)
                                source_completeness_2.append(src_completeness_2)
                                source_num_fibers.append(src_num_fibers)

                                source_sn.append(float(ri["sn"]))

                                if bid_shot == shot:
                                    miss_self.append(True)
                                else:
                                    miss_self.append(False)

                                #
                                # print("early break, check first miss")
                                # halt = True
                                # break

                            except:
                                print(f"Exception: {det} {bid_shot}", traceback.format_exc())
                                # halt = True
                                # break #(want to stop and take a look)

                            # valid_misses += 1

                    T_elixer_dets['num_matches'][sel_update] = mc_hits
                    T_elixer_dets['num_valid_misses'][sel_update] = valid_misses

                T_elixer_dets['num_valid_misses'][sel_update] = valid_misses


            #end for each det in a shot

            #intermediate write
            dtnow = datetime.now()
            dtnowstr = "%s:%s:%s.%s : " % (
            str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
            str(dtnow.microsecond).zfill(6))
            print(f"{dtnowstr} Intermediate write of cmplt table for group {group}")

            Cmplt_Table = Table(
                [source_group, source_shotid, source_detid, source_wave, source_line_flux, source_line_sigma,
                 source_num_fibers, source_sn, source_f50, source_completeness,
                 miss_shotid, miss_self, miss_found, miss_sn, miss_seeing, miss_tput, miss_num_fibers,
                 miss_f50, miss_completeness, ],
                names=['group', 'src_shotid', 'detectid', 'wave', 'line_flux', 'line_sigma',
                       'src_nfib', 'sn', 'src_f50', 'src_cmplt',
                       'bid_shotid', 'self', 'found', 'bid_sn', 'bid_seeing', 'bid_tput', 'bid_nfib',
                       'bid_f50', 'bid_cmplt'])

            Cmplt_Table.write(f"laurel_elixer_mc_completeness_{str(group_num).zfill(2)}.fits", overwrite=True,
                              format="fits")

        #end for each shot


Cmplt_Table = Table([source_group,source_shotid,source_detid,source_wave,source_line_flux,source_line_sigma,
                    source_num_fibers,source_sn, source_f50,source_completeness,
                    miss_shotid,miss_self,miss_found,miss_sn,miss_seeing,miss_tput,miss_num_fibers,
                    miss_f50,miss_completeness,],
                  names=['group','src_shotid','detectid','wave','line_flux','line_sigma',
                         'src_nfib','sn','src_f50','src_cmplt',
                         'bid_shotid','self','found','bid_sn','bid_seeing','bid_tput','bid_nfib',
                         'bid_f50','bid_cmplt'])

Cmplt_Table['wave'].format="%0.2f"
Cmplt_Table['line_flux'].format="%0.2f"
Cmplt_Table['line_sigma'].format="%0.2f"
Cmplt_Table['sn'].format="%0.2f"
Cmplt_Table['src_f50'].format="%0.2f"
Cmplt_Table['src_cmplt'].format="%0.4f"

Cmplt_Table['bid_sn'].format="%0.2f"
Cmplt_Table['bid_seeing'].format="%0.2f"
Cmplt_Table['bid_tput'].format="%0.4f"
Cmplt_Table['bid_f50'].format="%0.2f"
Cmplt_Table['bid_cmplt'].format="%0.4f"

dtnow = datetime.now()
dtnowstr = "%s:%s:%s.%s : " % (str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
                            str(dtnow.microsecond).zfill(6))
print(f"{dtnowstr} Building cmplt table for group {group}")

Cmplt_Table.write(f"laurel_elixer_mc_completeness_{str(group_num).zfill(2)}.fits",overwrite=True,format="fits")

dtnow = datetime.now()
dtnowstr = "%s:%s:%s.%s : " % (str(dtnow.hour).zfill(2), str(dtnow.minute).zfill(2), str(dtnow.second).zfill(2),
                            str(dtnow.microsecond).zfill(6))
print(f"{dtnowstr} Group {str(group_num).zfill(2)} done")

#then afterward, need to vstack the various tables