#!/usr/bin/env python3
# jwst_nircam_build_epsf_pooled.py
# Tested with: photutils 2.3.0

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import NDData
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

from photutils.psf import extract_stars, EPSFBuilder, EPSFStars
from photutils.background import Background2D, SExtractorBackground
from scipy.interpolate import RectBivariateSpline
from scipy.signal import convolve2d
import tqdm
from DavidsNM import miniNM_new, miniLM_new, save_img, save_patches
from FileRead import readcol
import glob
import sys

# ---------------------- user-tunable knobs ----------------------
CUTOUT_SIZE        = 21         # pixels; safe for NIRCam wings
OVERSAMPLING       = 2
MAX_ITERS          = 10
MIN_STARS_REQUIRED = 25         # minimum pooled cutouts to proceed
USE_DQ_MASK        = True
DQ_GOOD_VALUE      = 0          # mask anything != 0
KEEP_TOP_N_STARS   = 100        # cap on brightest pooled stars
BKG_BOX_SIZE       = None       # None => auto (≈ 1/20 of shorter image side, odd)
BKG_FILTER_SIZE    = (3, 3)
BKG_ESTIMATOR      = SExtractorBackground()


# ---------------------------------------------------------------

def save_patches_3d(threeD_item, name):
    save_patches([np.reshape(item, [len(item)*CUTOUT_SIZE, CUTOUT_SIZE]) for item in threeD_item], name)


def get_10x_design_matrix():
    x1d = np.arange(OVERSAMPLING*CUTOUT_SIZE)/float(OVERSAMPLING)
    x1d -= np.mean(x1d)

    x1d_10x = np.arange(CUTOUT_SIZE*10, dtype=np.float64)/10.
    x1d_10x -= np.mean(x1d_10x)

    design_matrix = np.zeros([(10*CUTOUT_SIZE)**2, (OVERSAMPLING*CUTOUT_SIZE)**2], dtype=np.float64)
    kernel = np.ones((10, 10), dtype=np.float64)

    for i in tqdm.trange((OVERSAMPLING*CUTOUT_SIZE)**2):
        P = np.zeros((OVERSAMPLING*CUTOUT_SIZE)**2, dtype=np.float64)
        P[i] = 1.
        
        P2d = np.reshape(P, [OVERSAMPLING*CUTOUT_SIZE]*2)
    
        ifn = RectBivariateSpline(x1d, x1d, P2d, kx = 2, ky = 2)
        result_2d = ifn(x1d_10x, x1d_10x)
        result_2d = convolve2d(result_2d, kernel, mode="same", boundary="symm")

        result_2d[0, :] = 0.
        result_2d[-1, :] = 0.
        result_2d[:, 0] = 0.
        result_2d[:, -1] = 0.
        
        design_matrix[:, i] = np.reshape(result_2d, (CUTOUT_SIZE*10)**2)
    return design_matrix

print("Getting design matrix")
x1d_10x = np.arange(CUTOUT_SIZE*10, dtype=np.float64)/10.
x1d_10x -= np.mean(x1d_10x)
design_matrix = get_10x_design_matrix()
#save_img(design_matrix, "design_matrix.fits")

def modelfn_all_stars(P, passdata):
    P2d = np.reshape(np.dot(design_matrix, P), [CUTOUT_SIZE*10]*2)
    #d_P2d = np.reshape(np.dot(design_matrix, d_P), [CUTOUT_SIZE*10]*2)
    
    ifn = RectBivariateSpline(x1d_10x, x1d_10x, P2d, kx = 2, ky = 2)
    #d_ifn = RectBivariateSpline(x1d_10x, x1d_10x, d_P2d, kx = 2, ky = 2)

    native_x = np.arange(CUTOUT_SIZE, dtype=np.float64)
    native_x -= np.mean(native_x)

    native_mod = np.zeros(passdata["data"].shape, dtype=np.float64)
    
    for i in range(len(passdata["data"])):
        native_mod[i] = (passdata["A"][i])*ifn(native_x - passdata["x0"][i], native_x - passdata["y0"][i]) + passdata["sky"][i]
    return native_mod


def residfn_all_stars(P, passdata):
    passdata = passdata[0]

    native_mod = modelfn_all_stars(P, passdata)
    resid = passdata["data"] - native_mod
    resid = resid.flatten()
    resid = resid[np.where(1 - np.isnan(resid))]
    
    
    return np.sign(resid) * np.sqrt(np.abs(resid))


def chi2fn_one_star(P, passdata):
    dat, ifn = passdata[0]
    native_x = np.arange(CUTOUT_SIZE, dtype=np.float64)
    native_x -= np.mean(native_x)

    mod = P[0]*ifn(native_x - P[1], native_x - P[2]) + P[3]
    resid = dat - mod
    resid = resid[np.where(1 - np.isnan(resid))]

    if np.abs(P[1]) > 5 or np.abs(P[2]) > 5:
        return 1e100
    
    return np.sum(np.abs(resid))

def fit_ampl_x0y0sky(passdata, ifn):
    for i in tqdm.trange(len(passdata["data"])):
        P, NA, NA = miniNM_new(ministart = [passdata["A"][i], passdata["x0"][i], passdata["y0"][i], passdata["sky"][i]],
                               miniscale = [passdata["A"][i]/10., 0.0, 0.0, 0.0],
                               chi2fn = chi2fn_one_star,
                               passdata = [passdata["data"][i], ifn], compute_Cmat = False)

        P, NA, NA = miniNM_new(ministart = [P[0], passdata["x0"][i], passdata["y0"][i], passdata["sky"][i]],
                               miniscale = [P[0]/10., 1.0, 1.0, 0.0],
                               chi2fn = chi2fn_one_star,
                               passdata = [passdata["data"][i], ifn], compute_Cmat = False)

        
        passdata["A"][i] = P[0]
        passdata["x0"][i] = P[1]
        passdata["y0"][i] = P[2]
        assert np.isclose(passdata["sky"][i], P[3])

        P, NA, NA = miniNM_new(ministart = [passdata["A"][i], passdata["x0"][i], passdata["y0"][i], passdata["sky"][i]],
                               miniscale = [0.0, 0.0, 0.0, 1.0],
                               chi2fn = chi2fn_one_star,
                               passdata = [passdata["data"][i], ifn], compute_Cmat = False)
        passdata["sky"][i] = P[3]

    passdata["sky"] -= np.median(passdata["sky"])
    return passdata

def compare_A(PSF1, PSF2):
    rel_A = np.zeros([10]*2, dtype=np.float64)
    for i in range(10):
        for j in range(10):
            rel_A[i,j] = np.sum(PSF1[i::10, j::10]*PSF2[i::10, j::10])
    rel_A /= np.median(rel_A)
    RMS = np.std(rel_A)
    print("rel_A", rel_A)
    print("RMS comparison", RMS)
    return RMS

def renorm_P_PSF(P_PSF):
    P2d = np.reshape(np.dot(design_matrix, P_PSF), [CUTOUT_SIZE*10]*2)
    norm_term = P2d[4::10, 4::10].sum()

    P2d /= norm_term
    P_PSF /= norm_term
    ifn = RectBivariateSpline(x1d_10x, x1d_10x, P2d, kx = 2, ky = 2)
    
    return P_PSF, P2d, ifn


    
def estimate_PSF_all_stars(all_dat_list, t_midpoints, PSF_rough_guess, verbose = False):
    
    passdata = dict(data = all_dat_list,
                    A = np.zeros(len(all_dat_list)) + 10000.,
                    x0 = np.zeros(len(all_dat_list)),
                    y0 = np.zeros(len(all_dat_list)),
                    sky = np.zeros(len(all_dat_list)))

    
    x1d = np.arange(OVERSAMPLING*CUTOUT_SIZE)/float(OVERSAMPLING)
    x1d -= np.mean(x1d)

    X, Y = np.meshgrid(x1d, x1d)

    P_PSF = np.zeros([OVERSAMPLING*CUTOUT_SIZE]*2, dtype=np.float64)
    
    for i in range(len(OVERSAMPLING)):
        for j in range(len(OVERSAMPLING)):
            P_PSF[i::OVERSAMPLING, j::OVERSAMPLING] = PSF_rough_guess
        
    P_PSF = np.reshape(P_PSF, (OVERSAMPLING*CUTOUT_SIZE)**2)
    P_PSF, P2d, ifn = renorm_P_PSF(P_PSF)
    
    
    passdata = fit_ampl_x0y0sky(passdata, ifn)
    the_mod = modelfn_all_stars(P = P_PSF, passdata = passdata)
    
    PSF_iter = 0

    if verbose:
        save_img(P2d, "PSF_10x_iter=%02i.fits" % 0)
        save_patches_3d(the_mod, "the_mod_iter=%02i.fits" % 0)
        save_patches_3d(all_dat_list - the_mod, "the_resid_iter=%02i.fits" % 0)

    
    
    last_PSF = P2d*1.

    
    while ((compare_A(P2d, last_PSF) > 0.005) or (PSF_iter < 2)) and (PSF_iter < 5):
        PSF_iter += 1
        last_PSF = P2d*1.
        
        passdata = fit_ampl_x0y0sky(passdata, ifn)
        
        P_PSF, NA, NA = miniLM_new(ministart = P_PSF,
                                   miniscale = np.ones((OVERSAMPLING*CUTOUT_SIZE)**2, dtype=np.float64),
                                   passdata = passdata,
                                   residfn = residfn_all_stars, verbose = True, maxiter = 10)

        P_PSF, P2d, ifn = renorm_P_PSF(P_PSF)
        
        passdata = fit_ampl_x0y0sky(passdata, ifn)

        the_mod = modelfn_all_stars(P = P_PSF, passdata = passdata)
        if verbose:
            save_img(P2d, "PSF_10x_iter=%02i.fits" % PSF_iter)
            save_patches_3d(the_mod, "the_mod_iter=%02i.fits" % PSF_iter)
            the_resid = all_dat_list - the_mod
            save_patches_3d(the_resid, "the_resid_iter=%02i.fits" % PSF_iter)

            the_norm_resid = the_resid*1.
            for i in range(len(the_resid)):
                the_norm_resid[i] /= the_mod[i].max()
            save_patches_3d(the_norm_resid, "the_norm_resid_iter=%02i.fits" % PSF_iter)

        
    return P_PSF, P2d

def _sky_to_pixel(wcs, ra_deg, dec_deg):
    sky = SkyCoord(ra_deg, dec_deg, unit='deg', frame='icrs')
    x, y = wcs.world_to_pixel(sky)
    return np.array(x), np.array(y)

def _in_bounds(x, y, nx, ny, margin):
    return (x > margin) & (x < nx - margin - 1) & (y > margin) & (y < ny - margin - 1)

def _background_2d(data, mask=None):
    ny, nx = data.shape
    if BKG_BOX_SIZE is None:
        box = max(32, int(min(nx, ny) // 20))
        if box % 2 == 0:
            box += 1
        box_size = (box, box)
    else:
        box_size = BKG_BOX_SIZE
    bkg = Background2D(
        data,
        box_size=box_size,
        filter_size=BKG_FILTER_SIZE,
        mask=mask,
        bkg_estimator=BKG_ESTIMATOR
    )
    return bkg.background



def _pixel_scale_arcsec(w: WCS):
    try:
        m = w.pixel_scale_matrix
        pixscale = np.sqrt((m**2).sum(axis=0)).mean() * 3600.0
        return float(pixscale)
    except Exception:
        return np.nan



def filter_close_pairs(x, y):
    # x, y: 1D arrays of star positions, same length
    x = np.asarray(x)
    y = np.asarray(y)

    coords = np.column_stack((x, y))
    n = len(coords)
    keep = np.ones(n, dtype=bool)

    min_sep = 15.0
    min_sep2 = min_sep**2  # compare squared distances to avoid sqrt

    for i in range(n):
        if not keep[i]:
            continue
        # Vector from star i to all later stars
        dx = coords[i+1:, 0] - coords[i, 0]
        dy = coords[i+1:, 1] - coords[i, 1]
        dist2 = dx*dx + dy*dy

        # Mark later stars within min_sep as "do not keep"
        close = dist2 < min_sep2
        keep[i+1:][close] = False

    # Filtered coordinates
    x_clean = x[keep]
    y_clean = y[keep]
    return x_clean, y_clean


def build_pooled_epsf(
        image_paths,
        ra_deg_array,
        dec_deg_array,
        cutout_size=CUTOUT_SIZE,
        oversampling=OVERSAMPLING,
        maxiters=MAX_ITERS,
        min_stars=MIN_STARS_REQUIRED,
        keep_top_n=KEEP_TOP_N_STARS,
        verbose = False):
    """Pool star cutouts across all images and build a single oversampled ePSF."""
    if len(image_paths) == 0:
        raise ValueError("Provide at least one _cal.fits image path.")

    all_dat_list = []
    pixscales = []

    for path in image_paths:
        print("Working on ", path)
        this_image_dat_list = []

        with fits.open(path) as hdul:
            print(hdul.info())
            
            sci = hdul['SCI'].data.astype(float)
            sat_mask = (hdul["GROUPDQ"].data & 2) != 0    # shape: (nint, ngroup, ny, nx)
            sci[sat_mask] = np.nan

            shdr = hdul['SCI'].header
            assert np.isclose(hdul[0].header["TGROUP"], 21.474)

            t_midpoints_one_ramp = (np.arange(len(sci[0]), dtype=np.float64)*2 + 1.5)*10.737

            t_midpoints = np.array(list(t_midpoints_one_ramp)*len(sci))
            
            print("t_midpoints", t_midpoints)

            
            
        with fits.open(path.replace("_uncallin.fits", "_tweakreg.fits")) as hdul:
            shdr = hdul['SCI'].header
            w = WCS(shdr)
            filter_name = hdul[0].header["FILTER"]
            print("filter_name", filter_name)
            out_name = path.split("_")[-2] + "_" + filter_name

            
            

        ny, nx = sci[0,0].shape
        print("nx", nx, "ny", ny)
        

        x, y = _sky_to_pixel(w, ra_deg_array, dec_deg_array)
        x, y = filter_close_pairs(x, y)
        
        keep = _in_bounds(x, y, nx, ny, margin=cutout_size // 2 + 2)
        x, y = x[keep], y[keep]
        if len(x) == 0:
            continue
        print("HERE4")

        print("x", x, len(x), "y", y)
        
        sci_diffs = sci[:, 1:] - sci[:, :-1]
        print("sci_diffs", sci_diffs.shape)
        sci_diffs = np.reshape(sci_diffs, [sci_diffs.shape[0]*sci_diffs.shape[1], sci_diffs.shape[2], sci_diffs.shape[3]])

        median_diff = np.nanmedian(sci_diffs, axis = 0)
        if verbose:
            save_img(sci_diffs, "sci_diffs.fits")
            save_img(median_diff, "median_diff.fits")

        
        
        for i in range(len(x)):
            this_image_dat_list.append(sci_diffs[:, int(np.around(y[i] - cutout_size/2)): int(np.around(y[i] + cutout_size/2)),
                                          int(np.around(x[i] - cutout_size/2)): int(np.around(x[i] + cutout_size/2))])

        for rough_iter in range(2):
            PSF_rough_guess = np.nanmedian(np.array([np.median(item, axis = 0) for item in this_image_dat_list]), axis = 0)
            PSF_rough_guess -= np.nanmedian(PSF_rough_guess)


            if verbose:
                save_img(PSF_rough_guess, "PSF_rough_guess_it=%i.fits" % rough_iter)

            flux_estimates = []
            for i in range(len(this_image_dat_list)):
                med_dat = np.nanmedian(this_image_dat_list[i], axis = 0)
                assert len(med_dat) == cutout_size

                med_dat -= np.nanmedian(med_dat)
                flux_estimates.append(np.nansum(med_dat*PSF_rough_guess))

            if rough_iter == 0:
                stars_to_keep = KEEP_TOP_N_STARS*2
            else:
                stars_to_keep = KEEP_TOP_N_STARS
                
            flux_cutoff = np.sort(flux_estimates)[-stars_to_keep]
            print("flux_cutoff", flux_cutoff)

            for i in range(len(this_image_dat_list))[::-1]:
                if flux_estimates[i] < flux_cutoff:
                    del this_image_dat_list[i]
            print("Remaining stars ", len(this_image_dat_list))

        if verbose:
            save_img(this_image_dat_list, "this_image_dat_list.fits")
            
        all_dat_list.extend(this_image_dat_list)
        pixscales.append(_pixel_scale_arcsec(w))
    print("HERE2")


    if len(all_dat_list) == 0:
        return {
            "status": "failed",
            "reason": "No usable star cutouts found across all images.",
            "images": image_paths
        }

    pixscales = np.array(pixscales)
    if len(pixscales) > 0 and (np.nanmax(pixscales) - np.nanmin(pixscales)) > 0.01:
        print("[WARN] Detected >0.01\" variation in pixel scale; avoid mixing SW/LW when pooling.")

    print("HERE1")

    if len(all_dat_list) < min_stars:
        return {
            "status": "failed",
            "reason": f"Too few pooled star cutouts ({len(pooled_list)} < {min_stars}).",
            "images": image_paths
        }

    print("HERE")


    all_dat_list = np.array(all_dat_list)
    
    if verbose:
        save_patches([np.nanmedian(item, axis = 0) for item in all_dat_list], "all_dat_list.fits")
        save_patches_3d(all_dat_list, "all_dat_list.fits")

    P, P2d = estimate_PSF_all_stars(all_dat_list, t_midpoints = t_midpoints, PSF_rough_guess = PSF_rough_guess, verbose = verbose)

    save_img(P2d, "PSF_10x_" + out_name.replace(".fits", "") + ".fits")


if __name__ == "__main__":
    # x_pix y_pix RA_deg Dec_deg
    [x_pix, y_pix, ra_deg, dec_deg] = readcol(sys.argv[1], 'ff,ff') # E.g., WD_jw02729001003_02105_nrca3.txt
    assert len(ra_deg) > 10

    verbose = sys.argv[2]
    images = sys.argv[3:]

    for image in images:
        assert image.count("uncallin.fits") == 1
        

    """
    for short_long in [0, 1]:
        glob_str = sys.argv[1].replace(".txt", "_cal.fits").replace("WD_", "")
        if short_long:
            the_glob = glob_str.replace("_nrc", "_*_nrc")
        else:
            assert glob_str.count("_nrca") or glob_str.count("_nrcb"), glob_str
            a_not_b = glob_str.count("_nrca")
            the_glob = glob_str.split("_nrc")[0] + "_*_nrc" + "a"*a_not_b + "b"*(1 - a_not_b) + "long_cal.fits"
            
        print("Glob ", the_glob)
        images = glob.glob(the_glob)
    """
    diag = build_pooled_epsf(
        images, ra_deg, dec_deg,
        oversampling=4,
        cutout_size=21,
        verbose = verbose,
    )
    
