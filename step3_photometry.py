# JWST NIRCam Level-3 detection + forced photometry (F150W2/F322W2)
# Goals: positions + rough colors (good enough to pick the WD cooling sequence)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, vstack, hstack
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
import astropy.units as u

from photutils.background import Background2D, SExtractorBackground
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

import sys

# --------------------------
# User inputs (edit paths)
# --------------------------
F150W2_SCI = sys.argv[1]
F150W2_ERR = sys.argv[2]

F322W2_SCI = sys.argv[3]
F322W2_ERR = sys.argv[4]

OUT_ECSV   = sys.argv[5]

# Detection tuning (conservative; tweak for depth vs. shredding)
NSIG_DET = 1.35          # threshold in sigma for segmentation detection
MIN_PIX  = 10            # min connected pixels for a detection
DEBLEND_NTHRESH = 64
DEBLEND_CONT    = 0.0005 # more aggressive deblending in crowded fields

# Background mesh (keep small-ish for crowding)
BKG_BOXSIZE   = 48
BKG_FILTERSIZE = 3

# Photometry: choose small apertures for S/N, then aperture-correct empirically
R_SMALL_FACTOR = 0.70    # r_small = 0.70 * FWHM (in pixels)
R_LARGE_SCALE  = 3.0     # r_large = 3.0 * r_small (for aperture correction)
ANN_IN_SCALE   = 4.0     # r_in = 4.0 * r_small
ANN_OUT_SCALE  = 6.0     # r_out = 6.0 * r_small
APCORR_TOP_N   = 300     # use N brightest *compact* detections to estimate apcorr

# Approximate PSF FWHM (diffraction), in arcsec, at wavelength lam [micron]
def jwst_fwhm_arcsec(lam_um):
    D = 6.5  # m
    return 1.03 * (lam_um * 1e-6 / D) * (206265.0)  # ~1.03 * λ/D in arcsec

# Pixel scale (arcsec/pix) from WCS
def pixscale_arcsec(w):
    # proj_plane_pixel_scales returns deg/pix on each axis
    from astropy.wcs.utils import proj_plane_pixel_scales
    sc = proj_plane_pixel_scales(w) * 3600.0  # arcsec/pix
    return float(np.mean(sc))

# Steradian per pixel (for flux conversion MJy/sr * sr -> MJy)
def pixarea_sr(hdr, w):
    if 'PIXAR_SR' in hdr:
        return float(hdr['PIXAR_SR'])
    # fallback from WCS scales (approx)
    sc = pixscale_arcsec(w) / 3600.0 * np.pi/180.0  # rad/pix along both axes (assuming square)
    return sc * sc

# AB mag from flux in Jy
def jy_to_abmag(f_jy):
    good = f_jy > 0
    m = np.full_like(f_jy, np.nan, dtype=float)
    m[good] = -2.5 * np.log10(f_jy[good]) + 8.90
    return m

# Build background, segmentation, and catalog on a given mosaic
def detect_on_mosaic(sci, wcs, nsig=1.3, minpix=10):
    bkg = Background2D(sci, box_size=BKG_BOXSIZE, filter_size=BKG_FILTERSIZE,
                       bkg_estimator=SExtractorBackground(), exclude_percentile=10.0)
    data_sub = sci - bkg.background
    thresh = bkg.background + nsig * bkg.background_rms

    segm = detect_sources(data_sub, thresh, npixels=minpix)
    if segm is None:
        raise RuntimeError("No detections found; lower threshold or check inputs.")
    segm = deblend_sources(data_sub, segm, npixels=minpix,
                           nlevels=DEBLEND_NTHRESH, contrast=DEBLEND_CONT)

    cat = SourceCatalog(data_sub, segm, wcs=wcs)
    tbl = cat.to_table()
    # Keep only compact-ish sources (very loose)
    if 'ellipticity' in tbl.colnames:
        tbl = tbl[tbl['ellipticity'] < 0.6]
    return tbl, bkg, data_sub, segm

# Aperture photometry (small + large for apcorr), with local annulus sky and err map
def forced_phot(table_master, sci, err, wcs, lam_um, band_label):
    # radii tied to FWHM (from λ/D) and pixel scale
    pscale = pixscale_arcsec(wcs)                   # arcsec/pix
    fwhm_as = jwst_fwhm_arcsec(lam_um)             # arcsec
    fwhm_pix = fwhm_as / pscale

    r_small = R_SMALL_FACTOR * fwhm_pix
    r_large = R_LARGE_SCALE  * r_small
    r_in    = ANN_IN_SCALE   * r_small
    r_out   = ANN_OUT_SCALE  * r_small

    positions = np.vstack((table_master['xcentroid'], table_master['ycentroid'])).T
    aper_s = CircularAperture(positions, r=r_small)
    aper_L = CircularAperture(positions, r=r_large)
    ann    = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

    # aperture areas in pixels
    n_pix_small = aper_s.area
    n_pix_large = aper_L.area

    # Local background (median in annulus)
    ann_masks = ann.to_mask(method='center')
    bkg_med = np.empty(len(positions), dtype=float)
    bkg_std = np.empty(len(positions), dtype=float)
    for i, m in enumerate(ann_masks):
        vals = m.multiply(sci)
        try:
            data = vals[m.data > 0]
        except:
            data = np.ones(2)
            
        if data.size < 10:
            bkg_med[i] = 0.0
            bkg_std[i] = np.nan
        else:
            _, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
            bkg_med[i] = median
            bkg_std[i] = std

    # Photometry in MJy/sr
    phot_s = aperture_photometry(sci, aper_s)['aperture_sum'].value
    phot_L = aperture_photometry(sci, aper_L)['aperture_sum'].value

    # Subtract local sky
    flux_s = phot_s - bkg_med * n_pix_small
    flux_L = phot_L - bkg_med * n_pix_large

    # Uncertainties: propagate err map + background estimation (approx)
    err_s = np.sqrt(aperture_photometry(err**2, aper_s)['aperture_sum'].value +
                    (n_pix_small * bkg_std)**2)
    err_L = np.sqrt(aperture_photometry(err**2, aper_L)['aperture_sum'].value +
                    (n_pix_large * bkg_std)**2)

    # Empirical aperture correction from bright, compact sources
    # Pick the top APCORR_TOP_N by flux_L (safer) among reasonably round sources
    idx_sort = np.argsort(flux_L)[::-1]
    idx_use = idx_sort[:min(APCORR_TOP_N, len(idx_sort))]
    # robust median of small/large flux ratio
    good = (flux_s[idx_use] > 0) & (flux_L[idx_use] > 0)
    ratio = flux_s[idx_use][good] / flux_L[idx_use][good]
    apcorr_mag = -2.5 * np.log10(np.nanmedian(ratio)) if np.any(good) else 0.0
    apcorr_fac = 10**(-0.4 * apcorr_mag)  # multiply small-aperture flux by this

    # Apply apcorr to small-aperture flux (preferred S/N)
    flux_corr = flux_s * apcorr_fac
    err_corr  = err_s * apcorr_fac

    # Convert to Jy: sum(MJy/sr) * pixel_area_sr * 1e6
    hdr = fits.getheader(F150W2_SCI if band_label.lower()=="f150w2" else F322W2_SCI)
    pix_sr = pixarea_sr(hdr, wcs)
    flux_jy = flux_corr * pix_sr * 1e6
    err_jy  = err_corr  * pix_sr * 1e6

    mag_ab = jy_to_abmag(flux_jy)
    emag   = np.full_like(mag_ab, np.nan, dtype=float)
    ok = (flux_jy > 0) & (err_jy > 0)
    emag[ok] = (2.5/np.log(10)) * (err_jy[ok] / flux_jy[ok])

    out = Table()
    out[f'{band_label}_r_small_pix'] = np.full(len(flux_jy), r_small)
    out[f'{band_label}_apcorr_mag']  = np.full(len(flux_jy), apcorr_mag)
    out[f'{band_label}_flux_jy']     = flux_jy
    out[f'{band_label}_fluxerr_jy']  = err_jy
    out[f'{band_label}_mag_ab']      = mag_ab
    out[f'{band_label}_magerr']      = emag
    return out

# --------------------------
# Load mosaics
# --------------------------
sw_hdul = fits.open(F150W2_SCI)
lw_hdul = fits.open(F322W2_SCI)
sw_img  = sw_hdul['SCI'].data
lw_img  = lw_hdul['SCI'].data
sw_wcs  = WCS(sw_hdul['SCI'].header)
lw_wcs  = WCS(lw_hdul['SCI'].header)

sw_err = fits.getdata(F150W2_ERR)  # MJy/sr 1-sigma per pixel
lw_err = fits.getdata(F322W2_ERR)

# --------------------------
# 1) Detect on F150W2
# --------------------------
print("Detecting on F150W2...")
det_tbl, sw_bkg, sw_sub, sw_segm = detect_on_mosaic(sw_img, sw_wcs, nsig=NSIG_DET, minpix=MIN_PIX)

print("det_tbl", det_tbl.columns)

# Keep a compact "stellar-like" subset for centroids; very permissive cuts
cols_keep = ['label', 'xcentroid','ycentroid','sky_centroid','kron_flux','area','eccentricity']
det_tbl = det_tbl[cols_keep]
# RA/Dec from SkyCoord column
sky = det_tbl['sky_centroid']
det_tbl['ra']  = sky.ra.deg
det_tbl['dec'] = sky.dec.deg

# --------------------------
# 2) Forced photometry on SW (F150W2) and LW (F322W2)
# --------------------------
print("Photometering F150W2...")
sw_phot = forced_phot(det_tbl, sw_img, sw_err, sw_wcs, lam_um=1.50, band_label='F150W2')

print("sw_phot", sw_phot)

print("Photometering F322W2 (forced at SW positions)...")
# Transform SW pixel centroids to LW pixel frame (use RA/Dec for robustness)
coords = SkyCoord(ra=det_tbl['ra']*u.deg, dec=det_tbl['dec']*u.deg, frame='icrs')
x_lw, y_lw = lw_wcs.world_to_pixel(coords)

# Make a copy of the master table with LW pixel positions for the forced apertures
master = det_tbl.copy()
master['xcentroid'] = x_lw
master['ycentroid'] = y_lw

print("lw_img", lw_img.shape, lw_img)
print("lw_err", lw_err.shape, lw_err)
print("master", master)


lw_phot = forced_phot(master, lw_img, lw_err, lw_wcs, lam_um=3.22, band_label='F322W2')

# --------------------------
# 3) Assemble final catalog
# --------------------------
out = Table()
out['id']     = det_tbl['label']
out['x_sw']   = det_tbl['xcentroid']   # F150W2 pixel coords
out['y_sw']   = det_tbl['ycentroid']
out['x_lw']   = x_lw                   # F322W2 pixel coords
out['y_lw']   = y_lw
out['ra']     = det_tbl['ra']
out['dec']    = det_tbl['dec']

out = hstack([out, sw_phot, lw_phot], join_type='exact')

# Colors convenient for WD selection
out['color_f150w2_minus_f322w2'] = out['F150W2_mag_ab'] - out['F322W2_mag_ab']

# Very simple quality flags
out['flag_snr_f150w2'] = (out['F150W2_flux_jy'] / np.maximum(out['F150W2_fluxerr_jy'], 1e-99)) > 5.0
out['flag_snr_f322w2'] = (out['F322W2_flux_jy'] / np.maximum(out['F322W2_fluxerr_jy'], 1e-99)) > 3.0
#out['flag_goodshape']  = (det_tbl['area'] < (np.pi * (3.0 * jwst_fwhm_arcsec(1.50) / pixscale_arcsec(sw_wcs))**2)) & (det_tbl['eccentricity'] < 0.5)

out.write(OUT_ECSV, format='ascii.ecsv', overwrite=True)
print(f"Done. Wrote {len(out)} rows to {OUT_ECSV}")
