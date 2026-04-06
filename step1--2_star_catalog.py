import numpy as np
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.ndimage import maximum_filter

# --- Inputs you already have ---
# 1. FITS file with the WCS of the original image
img_fname = sys.argv[1] #  E.g., F150W_stacked_resample.fits

# 2. Correlation image as a numpy array (same shape as original)
hdu = fits.open("my_corr.fits")
corr = hdu[0].data
hdu.close()

# --- Load WCS from the original image ---
hdu = fits.open(img_fname)
hdr = hdu["SCI"].header
wcs = WCS(hdr)
hdu.close()

# --- Find peaks in the correlation image ---
threshold = 0.75        # your correlation cutoff
box_size = 5           # neighborhood size for "local maximum" (in pixels)

# maximum_filter gives, for each pixel, the max in a box around it
local_max = maximum_filter(corr, size=box_size)

# A pixel is a peak if:
#  - its value equals the local maximum in its neighborhood
#  - and it's above the threshold
#  - and it's finite
peaks_mask = (corr == local_max) & (corr > threshold) & np.isfinite(corr)

# Get pixel indices of peaks
y_pix, x_pix = np.nonzero(peaks_mask)   # note: numpy returns (row, col) = (y, x)

# --- Convert pixel positions to RA/Dec ---
# WCS.pixel_to_world expects x, y in 0-based pixel coordinates,
# which matches numpy indices, so we can pass x_pix, y_pix directly.
skycoords = wcs.pixel_to_world(x_pix, y_pix)

ra_deg  = skycoords.ra.deg
dec_deg = skycoords.dec.deg

# Also grab the correlation values at the peaks
corr_vals = corr[y_pix, x_pix]

# --- Build a table and write to file ---
tab = Table(
    [
        x_pix,            # x pixel index (0-based)
        y_pix,            # y pixel index (0-based)
        ra_deg,           # RA in degrees
        dec_deg,          # Dec in degrees
        corr_vals         # correlation value
    ],
    names=["x_pix", "y_pix", "RA_deg", "Dec_deg", "corr"]
)

# Save as an ECSV (nice, human-readable + round-trippable), or use FITS if you prefer
tab.write("WD_candidates.ecsv", format="ascii.ecsv", overwrite=True)
# tab.write("star_candidates.fits", overwrite=True)

f = open("ds9.reg", 'w')
f.write("""# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
""")

for i in range(len(ra_deg)):
    f.write('circle(%f,%f,0.5")\n' % (ra_deg[i], dec_deg[i]))
f.close()

