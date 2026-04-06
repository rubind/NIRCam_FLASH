import numpy as np
from astropy.io import fits
import sys
import tqdm
from DavidsNM import save_patches, save_img
from photutils.detection import find_peaks
from skimage.feature import match_template


f = fits.open(sys.argv[1])
dat = f["SCI"].data*1.
f.close()

#f = fits.open(sys.argv[2])
#dat2 = f["SCI"].data*1.
#f.close()

NMAD = 1.4826*np.nanmedian(np.abs(dat - np.nanmedian(dat)))
print("NMAD", NMAD)

#NMAD2 = 1.4826*np.nanmedian(np.abs(dat2 - np.nanmedian(dat2)))
#print("NMAD2", NMAD2)

local_max = (
    (dat >= np.roll(dat, 1, 0)) &
    (dat >= np.roll(dat, -1, 0)) &
    (dat >= np.roll(dat, 1, 1)) &
    (dat >= np.roll(dat, -1, 1)))

inds = np.where((dat > 30*NMAD + np.nanmedian(dat))*(local_max))


print("inds", len(inds[0]), len(inds[0])/float(dat.shape[0]*dat.shape[1]))

all_cutouts = []
all_cutouts2 = []


for i in tqdm.trange(len(inds[0])):
    ind_i, ind_j = inds[0][i], inds[1][i]

    cutout = dat[ind_i - 7:ind_i+8,
                 ind_j - 7:ind_j+8]

    
    
    
    if np.any(np.isnan(cutout)):
        pass
    else:
        if np.isclose(cutout.max(), dat[ind_i, ind_j]):
            all_cutouts.append(cutout/cutout.max())
            #cutout = dat2[ind_i - 7:ind_i+8,
            #              ind_j - 7:ind_j+8]
            #all_cutouts2.append(cutout/cutout.max())

all_cutouts = np.array(all_cutouts)
save_patches(all_cutouts, "all_cutouts.fits")
est_PSF = np.median(all_cutouts, axis = 0)
save_img(est_PSF, "est_PSF.fits")

dat[np.where(np.isnan(dat) + np.isinf(dat))] = 0.

dat = dat

#corr = match_template(dat, est_PSF, pad_input=True) # DOESN'T WORK FOR SOME REASON!
#save_img(corr, "corr.fits")

my_corr = dat*0.

flat_PSF = est_PSF.flatten()



for i in tqdm.trange(7, len(dat) - 7):
    for j in range(7, len(dat[0]) - 7):
        cutout = dat[i - 7: i + 8,
                     j - 7: j + 8].flatten()
        my_corr[i,j] = np.corrcoef(flat_PSF, cutout)[0, 1]
save_img(my_corr, "my_corr.fits")

