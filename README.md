Step 0: Download _uncal.fits and _cal.fits files.

python ~/NIRCam_ramp/step1_resamp_all_ims.py (aligns images to Gaia, making _tweakreg.fits files then stacks each band, running more than one of the jobs at the same time may cause reference-file conflicts)

python ~/NIRCam_ramp/step1--2_find_stars.py F150W_stacked_resample.fits (computes approximate PSF, then correlation coefficient with that PSF over the image, no I don't know why the built in function for this has problems)
